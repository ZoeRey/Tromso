;  
;  Calculate energy transport for each Fourier wave number
;
;---------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;------------------------------------------------------------

begin

Nwave = 20

g=9.80665
Cp=1005.
p0= 100000.
PI = 3.14159
R=6371000.
L=2500000.

;----------------------------------------------------------
;Read data
;------------------------------------------------------------

path="/global/work/runegg/CESM/archive/SOM_ETout_CAM5/atm/hist/"
f0=path+"/SOM_ETout_CAM5.cam.h0.0087-07.nc/"
Tf0=addfile(f0,"r")
pathsave="../ARCTICLE/BACK-UP/"


Nlat=dimsizes(Tf0&lat)
Nlon=dimsizes(Tf0&lon)
Ntime=12


vEsetot = new((/Nwave+1,Ntime,Nlat/),float)
vQsetot = new((/Nwave+1,Ntime,Nlat/),float)

vEsetot!0="WaveNumb"
vEsetot!1="time"
vEsetot!2="lat"
vEsetot&WaveNumb = fspan(0,Nwave,Nwave+1)
vEsetot&time = fspan(0,11,12)
vEsetot&lat = Tf0&lat

vEsetot=0.
vQsetot=0.
copy_VarCoords(vEsetot,vQsetot)

delete(Tf0)

;----------------------------------------------------
;LOOP
;---------------------------------------------------

Nyear=1
Nmonths=12

do iyear=0, Nyear-1
  do imonth=0, Nmonths-1

if (imonth+1 .le. 9) then
 imx=sprinti("%0.2i",imonth+1)
else
 imx=tostring(imonth+1)
end if

    f1=path+"SOM_ETout_CAM5.cam.h0.008"+iyear+"-"+imx+".nc"
    Tf1=addfile(f1,"r")

    ai = Tf1->hyai
    bi = Tf1->hybi
    am = Tf1->hyam
    bm = Tf1->hybm
    ps = Tf1->PS
    U = Tf1->U
    T = Tf1->T
    V = Tf1->V
    Q = Tf1->Q
    Z = Tf1->Z3
  
    Nlev=dimsizes(am)
    lon=Tf1&lon
    lat=doubletofloat(Tf1&lat)
    coslat=cos(lat*2*PI/360.)
    d=2*PI*R*coslat 

;---------------------------------------------

    dP=dpres_hybrid_ccm(ps,p0,ai,bi)
    
    VdP = V*dP
    E=Cp*T+g*Z+0.5*(U*U+V*V)
 
    copy_VarCoords(VdP, E)

    Vtm = VdP(0,:,:,:)
    Etm = E(0,:,:,:)
    Qtm = Q(0,:,:,:)


;-----------------------------------------------
;Fourier Coefficients
;-----------------------------------------------

    av = new((/Nwave+1,Nlev,Nlat/),float)
    bv = new((/Nwave+1,Nlev,Nlat/),float)
    aE = new((/Nwave+1,Nlev,Nlat/),float)
    bE = new((/Nwave+1,Nlev,Nlat/),float)
    aQ = new((/Nwave+1,Nlev,Nlat/),float)
    bQ = new((/Nwave+1,Nlev,Nlat/),float)

    av!0="WaveNumb"
    av!1="lev"
    av!2="lat"
    av&WaveNumb = vEsetot&WaveNumb
    av&lev = V&lev
    av&lat = V&lat
    copy_VarCoords(av, bv)
    copy_VarCoords(av, aE)
    copy_VarCoords(av, bE)
    copy_VarCoords(av, aQ)
    copy_VarCoords(av, bE)


       do iwave=0,Nwave

          cos_Fourier = onedtond(doubletofloat(cos(2*PI*iwave*lon/360)),(/Nlev,Nlat,Nlon/))
          sin_Fourier = onedtond(doubletofloat(sin(2*PI*iwave*lon/360)),(/Nlev,Nlat,Nlon/))
 
          av(iwave,:,:) = dim_avg(Vtm*cos_Fourier)
          aE(iwave,:,:) = dim_avg(Etm*cos_Fourier)
          aQ(iwave,:,:) = dim_avg(Qtm*cos_Fourier)
          bv(iwave,:,:) = dim_avg(Vtm*sin_Fourier)
          bE(iwave,:,:) = dim_avg(Etm*sin_Fourier)
          bQ(iwave,:,:) = dim_avg(Qtm*sin_Fourier)
 
          delete(cos_Fourier)
          delete(sin_Fourier)

       end do
    

;--------------------------------------------------------------
;Waves
;--------------------------------------------------------------

    vEtmse = 2*(av*aE + bv*bE)
    vQtmse = 2*(av*aQ + bv*bQ)

    vEtmse(0,:,:)=av(0,:,:)*aE(0,:,:)+bv(0,:,:)*bE(0,:,:)
    vQtmse(0,:,:)=av(0,:,:)*aQ(0,:,:)+bv(0,:,:)*bQ(0,:,:)

    ;vEtmse(0,:,:)=vEtmse(0,:,:)/2
    ;vQtmse(0,:,:)=vEtmse(0,:,:)/2

    copy_VarCoords(av, vEtmse)
    copy_VarCoords(av, vQtmse)

    vEsetot(:,imonth,:)=vEsetot(:,imonth,:)+dim_sum_Wrap(vEtmse(WaveNumb|:,lat|:,lev|:)) 
    vQsetot(:,imonth,:)=vQsetot(:,imonth,:)+dim_sum_Wrap(vQtmse(WaveNumb|:,lat|:,lev|:))

    delete(vEtmse)
    delete(vQtmse)
print(imx)
   end do
 
 end do

;---------------------------------------------------------------------------------

vEsetot=vEsetot/Nyear/g
vQsetot=vQsetot/Nyear/g*L


do ilat=0 ,Nlat-1
vEsetot(:,:,ilat)=vEsetot(:,:,ilat)*1e-15*2*PI*R*cos(lat(ilat)*2*PI/360.)
end do

vEsetot_year=dim_sum_n(vEsetot,1)/Nmonths
vEsetot_year!1="lat"
vEsetot_year&lat=Tf1&lat

do ilat=0 ,Nlat-1
vQsetot(:,:,ilat)=vQsetot(:,:,ilat)*2*1e-15*PI*R*cos(lat(ilat)*2*PI/360.)
end do

vQsetot_year=dim_sum_n(vQsetot,1)/Nmonths
vQsetot_year!1="lat"
vQsetot_year&lat=Tf1&lat

Etot_Nwave=vEsetot_year+vQsetot_year
Etot_Nwave!1="lat"
Etot_Nwave&lat=Tf1&lat

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
delete(av)
delete(aE)
delete(aQ)
delete(bv)
delete(bE)
delete(bQ)


;==========================================
;Saving

system("/bin/rm -f "+pathsave+"vE_Nwave.nc")
system("/bin/rm -f "+pathsave+"vQ_Nwave.nc")
system("/bin/rm -f "+pathsave+"Etot_Nwave.nc")



ncdf1 = addfile(pathsave+"vE_Nwave.nc" ,"c")  
ncdf2 = addfile(pathsave+"vQ_Nwave.nc" ,"c")
ncdf3 = addfile(pathsave+"Etot_Nwave.nc","c")


ncdf1->vEsetot_year=vEsetot_year
ncdf2->vQsetot_year=vQsetot_year
ncdf3->Etot_Nwave=Etot_Nwave

print("DONE")

end
