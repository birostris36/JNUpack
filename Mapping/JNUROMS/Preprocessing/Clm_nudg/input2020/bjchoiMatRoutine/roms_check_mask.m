% grd = roms_get_grid(grd_file);

value = 1;

grdchk = grd;

grdchk.lon_psi(find(grd.mask_psi==value)) = [];
grdchk.lat_psi(find(grd.mask_psi==value)) = [];
grdchk.lon_rho(find(grd.mask_rho==value)) = [];
grdchk.lat_rho(find(grd.mask_rho==value)) = [];
grdchk.lon_u(find(grd.mask_u==value)) = [];
grdchk.lat_u(find(grd.mask_u==value)) = [];
grdchk.lon_v(find(grd.mask_v==value)) = [];
grdchk.lat_v(find(grd.mask_v==value)) = [];

plot(grdchk.lon_psi,grdchk.lat_psi,'ro',...
    grdchk.lon_rho,grdchk.lat_rho,'gx',...
    grdchk.lon_u,grdchk.lat_u,'b+',...
    grdchk.lon_v,grdchk.lat_v,'cd')
legend('psi','rho','u','v')
