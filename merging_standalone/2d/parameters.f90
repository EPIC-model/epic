module parameters

 !Number of x & z grid divisions:
integer,parameter:: nx=40, nz=20

 !Maximum number of parcels:
integer,parameter:: nm=64*nx*nz

 !Box width and height:
double precision,parameter:: ellx=4.d0, ellz=2.d0

end module parameters
