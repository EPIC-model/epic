module parameters

 !Number of x, y and z grid divisions:
integer,parameter:: nx=160, ny=80, nz=40

 !Maximum number of parcels:
integer,parameter:: nm=64*nx*ny*nz

 !Box width, breadth and height:
double precision,parameter:: ellx=4.d0, elly=2.d0, ellz=1.d0

end module parameters
