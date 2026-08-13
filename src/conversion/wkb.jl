# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -------------------------------------
# WKB (Well-Known Binary) -> Meshes.jl
# -------------------------------------

# supports standard, extended, and ISO WKB geometry with Z dimensions (M/ZM not supported)
function wkb2meshes(buff, crs)
  # swap bytes of coordinates if necessary
  swapbytes = isone(read(buff, UInt8)) ? ltoh : ntoh

  # retrieve WKB geometry type
  wkbtype = read(buff, UInt32)

  # SQL/MM Part 3 and SFSQL 1.2 use offsets to
  # indicate the presence of higher dimensional
  # coordinates in a WKB geometry.
  #
  # `hasz` reflects what THIS geometry record actually declares in its type
  # code, which is the authoritative source for the on-disk layout. It can
  # disagree with `crs` (e.g. a Projected/EPSG CRS built from
  # gpkg_geometry_columns metadata has no 3D variant): in that case we must
  # still consume the Z bytes below to stay aligned with the rest of the
  # record, even though the resulting point can't carry that Z value.
  hasz = false
  if wkbtype ≥ 1001 && wkbtype ≤ 1007
    # 1000 (Z)
    wkbtype -= UInt32(1000)
    hasz = true
  elseif wkbtype ≥ 2001 && wkbtype ≤ 2007
    # 2000 (M)
    wkbtype -= UInt32(2000)
  elseif wkbtype ≥ 3001 && wkbtype ≤ 3007
    # 3000 (ZM)
    wkbtype -= UInt32(3000)
    hasz = true
  elseif wkbtype > 0x80000000
    # 99-402 was a short-lived extension to SFSQL 1.1
    # that used a high-bit flag to indicate the presence
    # of Z coordinates in a WKB geometry
    wkbtype -= 0x80000000
    hasz = true
  elseif wkbtype > 0x40000000
    # The M coordinate value allows the application environment
    # to associate some measure with the point values
    # this high-bit flag indicates the presence of M dimension
    wkbtype -= 0x40000000
  end

  # convert WKB geometry type to Meshes.jl type
  if wkbtype == 1
    wkb2point(buff, crs, swapbytes, hasz)
  elseif wkbtype == 2
    wkb2chain(buff, crs, swapbytes, hasz)
  elseif wkbtype == 3
    wkb2poly(buff, crs, swapbytes, hasz)
  elseif 4 ≤ wkbtype ≤ 7
    # do a recursive call to read inner geometries
    ngeoms = read(buff, UInt32)
    geoms = [wkb2meshes(buff, crs) for _ in 1:ngeoms]
    Multi(geoms)
  else
    error("Unsupported WKB Geometry Type: $wkbtype")
  end
end

wkb2point(buff, crs, swapbytes, hasz) = Point(wkb2coords(buff, crs, swapbytes, hasz))

wkb2points(buff, npoints, crs, swapbytes, hasz) = [wkb2point(buff, crs, swapbytes, hasz) for _ in 1:npoints]

function wkb2chain(buff, crs, swapbytes, hasz)
  npoints = read(buff, UInt32)
  points = wkb2points(buff, npoints, crs, swapbytes, hasz)
  if first(points) == last(points)
    while first(points) == last(points) && length(points) ≥ 2
      pop!(points)
    end
    Ring(points)
  else
    Rope(points)
  end
end

function wkb2poly(buff, crs, swapbytes, hasz)
  nrings = read(buff, UInt32)
  rings = [wkb2chain(buff, crs, swapbytes, hasz) for _ in 1:nrings]
  PolyArea(rings)
end

# Number of coordinates the CRS itself can represent (2 for LatLon/Projected/
# Cartesian2D, 3 for LatLonAlt/Cartesian3D).
_crsncoords(crs) = CoordRefSystems.ncoords(crs)

function wkb2coords(buff, crs, swapbytes, hasz)
  n = _crsncoords(crs)
  xy = ntuple(min(n, 2)) do _
    swapbytes(read(buff, Float64))
  end
  z = if hasz
    swapbytes(read(buff, Float64))
  else
    nothing
  end

  if crs <: LatLon
    crs(xy[2], xy[1])
  elseif crs <: LatLonAlt
    crs(xy[2], xy[1], something(z, 0.0))
  elseif n == 3
    crs(xy[1], xy[2], something(z, 0.0))
  else
    crs(xy...)
  end
end

# -------------------------------------
# Meshes.jl -> WKB (Well-Known Binary)
# -------------------------------------

_wkbtype(::Point) = 0x00000001
_wkbtype(::Chain) = 0x00000002
_wkbtype(::Polygon) = 0x00000003
_wkbtype(::MultiPoint) = 0x00000004
_wkbtype(::MultiChain) = 0x00000005
_wkbtype(::MultiPolygon) = 0x00000006
_wkbtype(::Multi) = 0x00000007

function meshes2wkb!(buff, geom)
  wkbtype = _wkbtype(geom)

  # wkbByteOrder = Little Endian
  write(buff, one(UInt8))

  # wkbGeometryType
  ncoords = CoordRefSystems.ncoords(crs(geom))
  ncoords == 3 ? write(buff, wkbtype + UInt32(1000)) : write(buff, wkbtype)

  if 1 ≤ wkbtype ≤ 3
    _meshes2wkb!(buff, geom)
  elseif 4 ≤ wkbtype ≤ 7
    gs = parent(geom)
    write(buff, UInt32(length(gs)))
    for g in gs
      meshes2wkb!(buff, g)
    end
  else
    throw(ErrorException("Well-Known Binary Geometry unknown: $wkbtype"))
  end
end

_meshes2wkb!(buff, point::Point) = _meshes2wkb!(buff, coords(point))

function _meshes2wkb!(buff, c::LatLon)
  write(buff, htol(ustrip(c.lon)))
  write(buff, htol(ustrip(c.lat)))
end

function _meshes2wkb!(buff, c::LatLonAlt)
  write(buff, htol(ustrip(c.lon)))
  write(buff, htol(ustrip(c.lat)))
  write(buff, htol(ustrip(c.alt)))
end

function _meshes2wkb!(buff, c::Projected)
  write(buff, htol(ustrip(c.x)))
  write(buff, htol(ustrip(c.y)))
end

function _meshes2wkb!(buff, c::Cartesian2D)
  write(buff, htol(ustrip(c.x)))
  write(buff, htol(ustrip(c.y)))
end

function _meshes2wkb!(buff, c::Cartesian3D)
  write(buff, htol(ustrip(c.x)))
  write(buff, htol(ustrip(c.y)))
  write(buff, htol(ustrip(c.z)))
end

function _meshes2wkb!(buff, chain::Chain)
  npoints = nvertices(chain)
  points = vertices(chain)
  if isclosed(chain)
    write(buff, UInt32(npoints + 1))
    for point in points
      _meshes2wkb!(buff, point)
    end
    _meshes2wkb!(buff, first(points))
  else
    write(buff, UInt32(npoints))
    for point in points
      _meshes2wkb!(buff, point)
    end
  end
end

function _meshes2wkb!(buff, poly::Polygon)
  rs = rings(poly)
  write(buff, UInt32(length(rs)))
  for r in rs
    _meshes2wkb!(buff, r)
  end
end
