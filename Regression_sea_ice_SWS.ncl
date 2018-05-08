begin

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;-----------------------------------------------------------
;DEFINE VARIABLE


datapath_sea_ice="/global/work/runegg/interim/surface/2.5x2.5/an/SIC/"
datapath_SWS="/global/work/runegg/interim/surface/2.5x2.5/fc.daily.24h/SWS/"
Tf0=addfile(datapath_sea_ice+"SIC.2017.12.nc","r")
Tf01=addfile(datapath_SWS+"SWS.2017.12.nc","r")

;SIC=Tf0->CI_GDS0_SFC
;SWS=Tf01->SSR_GDS0_SFC
;printVarSummary(SIC)
;printVarSummary(SWS)

PI=3.14159
R=6371000.

Nyear=13
Nmonths=12
Nstat=200

Nlat=dimsizes(Tf0&g0_lat_1)
Nlon=dimsizes(Tf0&g0_lon_2)

lat=doubletofloat(Tf0&g0_lat_1)
lon=doubletofloat(Tf0&g0_lon_2)

Sum_cos_lat=0.
do ilat=0,Nlat-1
 Sum_cos_lat=Sum_cos_lat+cos(lat(ilat)*0.0174533)
end do


SIC=new((/365*Nyear,Nlat,Nlon/),float)
SIC!0="time"
SIC!1="lat"
SIC!2="lon"
SIC&time=fspan(0,364*Nyear,365*Nyear)
SIC&lat=Tf0&g0_lat_1
SIC&lon=Tf0&g0_lon_2
SIC=0.

SWS=new((/365*Nyear,Nlat,Nlon/),float)
copy_VarCoords(SIC,SWS)
SWS=0.

;-----------------------------------------------------------------
;READ SEA-ICE

month=(/31,28,31,30,31,30,31,31,30,31,30,31/)
monthb=(/0,31,59,90,120,151,181,212,243,273,304,334/)
monthe=(/30,58,89,119,150,180,211,242,272,303,333,364/)

do iyear=0,Nyear-1
do imonth=0, Nmonths-1

  if (imonth+1 .le. 9) then
   imx=sprinti("%0.2i",imonth+1)
  else
   imx=tostring(imonth+1)
  end if
;print(imx)

  if(iyear .le. 20) then
   iyx=tostring(1979+iyear)
  else
   iyx=tostring(2000+(iyear-21))
  end if

  datapath_sea_ice="/global/work/runegg/interim/surface/2.5x2.5/an/SIC/"
  Tf1=addfile(datapath_sea_ice+"SIC."+iyx+"."+imx+".grb","r")

  SIC_month=Tf1->CI_GDS0_SFC

;printMinMax(SIC_month,0)

  SIC(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:,:)=SIC_month(0:month(imonth)-1,:,:)

delete(SIC_month)

end do
end do

;print(SIC)

S=new((/Nlat,Nlon/),float)
S=0.
do ilat=0,Nlat-1
 do ilon=0, Nlon-1
      if (.not.ismissing(SIC(0,ilat,ilon)))
        S(ilat,ilon)=1
end if
end do
end do 


;-----------------------------------------------------------
;READ SOLAR RAIDATION


do iyear=0,Nyear-1
do imonth=0, Nmonths-1

  if (imonth+1 .le. 9) then
   imx=sprinti("%0.2i",imonth+1)
  else
   imx=tostring(imonth+1)
  end if

  if(iyear .le. 20) then
   iyx=tostring(1979+iyear)
  else
   iyx=tostring(2000+(iyear-21))
  end if

   Tf1=addfile(datapath_SWS+"SWS."+iyx+"."+imx+".grb","r")

   SWS_time=Tf1->SSR_GDS0_SFC

;printMinMax(SWS_time,0)


   SWS(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:,:)=SWS_time(0:month(imonth)-1,:,:)


delete(SWS_time)

end do
end do

SWS=SWS/86400.

;printVarSummary(SWS)


;-----------------------------------------------------------
;DAILY CLIMATOLOGY

yyyyddd=new((/365*Nyear/),integer)
ddd=ispan(1,365,1)

do iyear=0,Nyear-1
 year0=fspan(1979+iyear,1979+iyear,365)
 year=floattoint(year0)
 yyyyddd(iyear*365:(iyear+1)*365-1)=year*1000+ddd 
end do  

;printVarSummary(yyyyddd)

SWS_clim=clmDayTLL(SWS,yyyyddd)         
SIC_clim=clmDayTLL(SIC,yyyyddd)

SWS_sm=runave_n(SWS_clim,8,-1,0)
SIC_sm=runave_n(SIC_clim,8,-1,0)

SWS_an=calcDayAnomTLL(SWS,yyyyddd,SWS_sm)
SIC_an=calcDayAnomTLL(SIC,yyyyddd,SIC_sm)

SWS_an!0="time"
SWS_an&time=fspan(0,Nyear*365-1,365*Nyear)

;printVarSummary(SWS_an)
;printMinMax(vQp_an,0)
;printMinMax(SIC_an,0)  

delete(SWS_sm)
delete(SWS_clim)
delete(SIC_sm)
delete(SIC_clim)

;-----------------------------------------------------------
;REGRESSION FOR THE ENTIRE YEAR

lag=40
Nlag=2*lag+1

  SWS_reg0=dim_avg(SWS_an(:,:,:))
  SWS_reg=dim_sum(SWS_reg0(:,0:24))
;printVarSummary(SWS_reg)

  Reg=new((/Nlag,Nlat,Nlon/),float)
  Reg!0="time_lag"
  Reg!1="lat"
  Reg!2="lon"
  Reg&time_lag=fspan(-lag,lag,2*lag+1)
  Reg&lat=Tf0&g0_lat_1
  Reg&lon=Tf0&g0_lon_2
  Reg=0.

  N=dimsizes(SWS_an&time)

  do ilat=0,Nlat-1
   do ilon=0,Nlon-1
    do ilag=-lag,lag

  Reg(lag+ilag,ilat,ilon)=avg(SIC_an(lag+ilag:N-lag-1+ilag,ilat,ilon)*SWS_reg(lag:N-1-lag))/(variance(SWS_reg(lag:N-1-lag)))


    end do
   end do
  end do


;printMinMax(Reg,0)
;print(Reg&lat)
;print(Reg&lon)

;do ilag=0,Nlag-1
;do ilat=0,Nlat-1
;do ilon=0,Nlon-1
; Reg(ilag,ilat,ilon)=Reg(ilag,ilat,ilon)*R*R*cos(lat(ilat)*0.0174533)*(2.5*0.0174533)*(2.5*0.0174533)
;end do
;end do
;end do


Reg_plot0=dim_avg(Reg(:,0:24,:))     ;30 to 90 °N
;Reg_plot=dim_avg(Reg_plot0)
lat_arctic=(/90,87.5,85,82.5,80,77.5,75,72.5,70,67.5,65,62.5,60,57.5,55,52.5,50,47.5,45,42.5,40,37.5,35,32.5,30/)
 do ilat=0,24
Reg_plot0(:,ilat)=Reg_plot0(:,ilat)*cos(lat_arctic(ilat)*0.0174533)
 end do
Reg_plot=dim_sum(Reg_plot0)/Sum_cos_lat

Reg_plot!0="time_lag"
Reg_plot&time_lag=fspan(-lag,lag,2*lag+1)

printMinMax(Reg_plot,0)

delete(Reg_plot0)
;------------------------------------------------------------
;STATISTICAL SIGNIFICANCE  

;RANDOM PHASE FOR vQp_reg


Stat_sign=new((/Nlag/),float)
Stat_sign!0="time_lag"
Stat_sign&time_lag=fspan(-lag,lag,2*lag+1)
;Stat_sign&lat=Tf0&g0_lat_1
Stat_sign=0.

SWS_reg_rand_phase=new(Nyear*365,float)
SWS_reg_rand_phase!0="time"
SWS_reg_rand_phase&time=fspan(0,365*Nyear-1,Nyear*365)

;file_back = "../REGRESSION/STAT_SIGN/Stat_sign_sea_ice.nc"
;a=addfile("../REGRESSION/STAT_SIGN/Stat_sign.nc","r")
;a = addfile(file_back,"r")
;Stat_sign = a->$"Stat_sign"$ ; time,lat,lon
;iSS= a->$"iSS"$
iSS=0
;print(iSS)

 do i=iSS,Nstat-1


 x=generate_unique_indices(365*Nyear)
 do k=0,365*Nyear-1
  SWS_reg_rand_phase(k)=SWS_reg(x(k))
 end do

SWS_reg_rand_phase!0="time"
SWS_reg_rand_phase&time=fspan(0,365*Nyear-1,365*Nyear)

;RGRESSION WITH vQp_reg_rand_phase

  Reg_rand_phase=new((/Nlag,Nlat,Nlon/),float)
  copy_VarCoords(Reg,Reg_rand_phase)

  N=dimsizes(SWS_reg_rand_phase&time)

  do ilat=0,Nlat-1
   do ilon=0,Nlon-1
    do ilag=-lag,lag
  Reg_rand_phase(lag+ilag,ilat,ilon)=avg(SIC_an(lag+ilag:N-lag-1+ilag,ilat,ilon)*SWS_reg_rand_phase(lag:N-1-lag))/(variance(SWS_reg_rand_phase(lag:N-1-lag)))

    end do
   end do
  end do

printVarSummary(Reg_rand_phase)

;do ilag=0,Nlag-1
;do ilat=0,Nlat-1
; do ilon=0,Nlon-1
; Reg_rand_phase(ilag,ilat,ilon)=Reg_rand_phase(ilag,ilat,ilon)*R*R*cos(lat(ilat)*0.0174533)*(2.5*0.0174533)*(2.5*0.0174533)
; end do
;end do
;end do

 Reg_rand_phase2=dim_avg(Reg_rand_phase(:,0:24,:))     ;30 to 90 °N
 do ilat=0,24
Reg_rand_phase2(:,ilat)=Reg_rand_phase2(:,ilat)*cos(lat_arctic(ilat)*0.0174533)
 end do
Reg_rand_phase_plot=dim_sum(Reg_rand_phase2)/Sum_cos_lat


;COMPARISON

  Stat_sign=Stat_sign+where(abs(Reg_rand_phase_plot).lt.abs(Reg_plot),1.,0.)



;  system("rm -f "+file_back)
;  ncdfRaw = addfile(file_back ,"c")
;  fAtt = True
;  fAtt@title = "MC significance test"
;  fAtt@Conventions = "None"
;  fAtt@creation_date = systemfunc("date")
;  fileattdef( ncdfRaw, fAtt )
;  ncdfRaw->Stat_sign = Stat_sign
;  ncdfRaw->iSS = i+1



delete(Reg_rand_phase)
delete(Reg_rand_phase2)
print(i)
end do


Stat_sign=Stat_sign/Nstat*100


;---------------------------------------------------------
;PLOT


Stat_sign_plot=Stat_sign
Stat_sign_plot!0="time_lag"
Stat_sign_plot&time_lag=fspan(-lag,lag,2*lag+1)

printMinMax(Stat_sign_plot,0)

Area=0.
 do ilat=0,Nlat-1
  do ilon=0,Nlon-1
 Area=Area+(R*R*cos(lat(ilat)*0.0174533)*(2.5*0.0174533)*(2.5*0.0174533))*S(ilat,ilon)
end do
end do

Reg_split=new((/3,Nlag/),float)

Reg_split(0,:)=Reg_plot*Area
Reg_split(1,:)=where(Stat_sign_plot .gt. 95,Reg_plot*Area,Reg_plot@_FillValue)
;Reg_split(2,:)=where(Stat_sign_plot .gt. 99,Reg_plot*Area,Reg_plot@_FillValue)
Reg_split(2,:)=0.

;----------------------------------------------------------
;PLOT 2D


wks  = gsn_open_wks("eps","regression_sea_ice_SWS")
 

res = True 
res@tiMainString="SSR" 
res@tiXAxisString = "Time lag (days)"  
res@tiYAxisString   = "Regression coefficient (m^2/W)"  
res@xyLineColors=(/"black","green","black"/)
res@xyLineThicknesses=(/1,4,1/)
res@xyDashPatterns=(/0,0,1/)
     
plot=gsn_csm_xy(wks,Reg_plot&time_lag, Reg_split,res)   

;-------------------------------------------------------------
;SAVE

path_saving = "../SOLAR_RADIATION/"
system("/bin/rm -f "+path_saving+"Sea_ice_SWS.nc")
nc = addfile(path_saving + "Sea_ice_SWS.nc","c")
nc->Reg_split_SWS=Reg_split



;-----------------------------------------------------------------
print("REGRESSION SEA ICE DONE")

end

