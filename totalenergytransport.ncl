;  Calculate total energy transport
;
;---------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;------------------------------------------------------------

begin

Nwave = 20

PI = 3.14159
L = 2500000.
g=9.80665
Cp=1005.
p0= 100000.
R=6371000.

;----------------------------------------------------------
;Read data


path="/global/work/runegg/CESM/archive/SOM_ETout_CAM5/atm/hist/"
f0=path+"/SOM_ETout_CAM5.cam.h0.0087-07.nc/"
Tf0=addfile(f0,"r")
pathsave="../ARCTICLE/BACK-UP/"

Nlat=dimsizes(Tf0&lat)
Nlon=dimsizes(Tf0&lon)
Ntime=12


;Total

vEsum = new((/12,Nlat/),float)   ;create a float
vQsum = new((/12,Nlat/),float)

vEsum!0="time"                   ;name the first dimension
vEsum!1="lat"
vEsum&time = fspan(0,11,12)      ;declare the first dimension
vEsum&lat = Tf0&lat

vEsum=0.
vQsum=0.

copy_VarCoords(vEsum,vQsum)


delete(Tf0)
;----------------------------------------

Nyear=2
Nmonths=12

do iyear=0, Nyear-1            ;for each year, for each month ...
  do imonth=0, Nmonths-1

if (imonth+1 .le. 9) then
 imx=sprinti("%0.2i",imonth+1)
else
 imx=tostring(imonth+1)
end if
 
    f1=path+"SOM_ETout_CAM5.cam.h0.008"+iyear+"-"+imx+".nc"    
    Tf1=addfile(f1,"r")

    ai = Tf1->hyai      ;in the file, recover some variables
    bi = Tf1->hybi
    am = Tf1->hyam
    bm = Tf1->hybm
    ps = Tf1->PS
    U = Tf1->U
    T = Tf1->T
    V = Tf1->V
    Q = Tf1->Q
    Z = Tf1->Z3
   
    Ntime=dimsizes(Tf1&time)
    Nlev=dimsizes(am)
    lon=Tf1&lon
    lat=doubletofloat(Tf1&lat)
    
      
;-------------------------------------------------------

    dP=dpres_hybrid_ccm(ps,p0,ai,bi)
    
    VdP = V*dP
    E=Cp*T+g*Z+0.5*(U*U+V*V)

    copy_VarCoords(VdP,E)

    Vtm = VdP(0,:,:,:)
    Etm = E(0,:,:,:)
    Qtm = Q(0,:,:,:)

    Energy_DS = (/dim_avg_Wrap(Vtm*Etm/g)/) ;calculate dry-static energy
    Energy_L = (/dim_avg_Wrap(Qtm*Vtm*L/g)/) ;calculate latent energy

    Energy_DS!0 = "lev"
    Energy_DS!1 = "lat"
    Energy_DS&lev = Tf1&lev
    Energy_DS&lat = Tf1&lat
    Energy_L!0 = "lev"
    Energy_L!1 = "lat"

    vEsum(imonth,:) = vEsum(imonth,:)+dim_sum_Wrap(Energy_DS(lat|:,lev|:))  
;calculate the sum of the DS_Energy for each month of each year
    vQsum(imonth,:) = vQsum(imonth,:)+ dim_sum_Wrap(Energy_L(lat|:,lev|:))
    
        end do
 end do

;-------------------------------------------------------
printVarSummary(vEsum)

vEsum=vEsum/Nyear        ;divide the sum by the number of year
do ilat=0,Nlat-1
vEsum(:,ilat)=vEsum(:,ilat)*2*PI*R*cos(lat(ilat)*2*PI/360.)/1e15   
end do

vQsum=vQsum/Nyear
do ilat=0,Nlat-1
vQsum(:,ilat)=vQsum(:,ilat)*2*PI*R*cos(lat(ilat)*2*PI/360.)/1e15
end do

vE=dim_avg_n(vEsum,0)     ;do an average over a year
vE=vE(::-1)
vQ=dim_avg_n(vQsum,0)
vQ=vQ(::-1)

Etot=dim_sum_n(vEsum,0)+dim_sum_n(vQsum,0)
;calculate the total energy transport
Etot=Etot(::-1)
Etot!0="lat"
Etot&lat=Tf1&lat
Etot=Etot/Nmonths


;------------------------------------------------------------------
;DELETE

delete(U)
delete(T)
delete(V)
delete(Q)
delete(Z)
delete(dP)
delete(E)
delete(Vtm)
delete(Etm)
delete(Qtm)
delete(Energy_DS)
delete(Energy_L)

;---------------------------------------------------------------------
;SAVING

system("/bin/rm -f "+pathsave+"vE.nc")
system("/bin/rm -f "+pathsave+"vQ.nc")
system("/bin/rm -f "+pathsave+"Etot.nc")

ncdf1 = addfile(pathsave + "vE.nc","c")
ncdf2 = addfile(pathsave + "vQ.nc","c")
ncdf3 = addfile(pathsave + "Etot.nc","c")

ncdf1->vE=vE
ncdf2->vQ=vQ
ncdf3->Etot=Etot

print("DONE")




end
