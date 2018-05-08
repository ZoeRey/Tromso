begin

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;-----------------------------------------------------------
;DEFINE VARIABLE


datapath_sea_ice="/global/work/runegg/interim/surface/2.5x2.5/an/SIC/"
datapath_energy="/global/work/runegg/interim/EnergySplit/Waves/"
Tf0=addfile(datapath_sea_ice+"SIC.2017.12.nc","r")
Tf01=addfile(datapath_energy+"vQtot.1980.01.WN20.nc","r")

PI=3.14159
R=6371000.

Nyear=33
Nmonths=12
Nstat=200

Nlat0=dimsizes(Tf0&g0_lat_1)
Nlat01=dimsizes(Tf01&lat)
Nlon=dimsizes(Tf0&g0_lon_2)
Nwave=dimsizes(Tf01&WaveNumb)
lat=doubletofloat(Tf01&lat)
lat0=doubletofloat(Tf0&g0_lat_1)
lon=doubletofloat(Tf0&g0_lon_2)

Sum_cos_lat=0.
do ilat=0,Nlat0-1
 Sum_cos_lat=Sum_cos_lat+cos(lat(ilat)*0.0174533)
end do

SIC=new((/365*Nyear,Nlat0,Nlon/),float)
SIC!0="time"
SIC!1="lat"
SIC!2="lon"
SIC&time=fspan(0,364*Nyear,365*Nyear)
SIC&lat=Tf0&g0_lat_1
SIC&lon=Tf0&g0_lon_2
SIC=0.

vQ=new((/365*Nyear,Nlat01/),float)
vQ!0="time"
vQ!1="lat"
vQ&time=fspan(0,364*Nyear,365*Nyear)
vQ&lat=Tf01&lat
vQ=0.

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
  datapath2_sea_ice="../nc/"

  fileogrb = datapath_sea_ice+"SIC."+iyx+"."+imx+".grb"
  fileonc = datapath2_sea_ice+"SIC."+iyx+"."+imx+".nc"
  system("cdo -R -f nc -t ecmwf copy "+fileogrb+" "+fileonc)
  system("/bin/rm -f "+fileonc)

  Tf1=addfile(datapath_sea_ice+"SIC."+iyx+"."+imx+".grb","r")

  SIC_month=Tf1->CI_GDS0_SFC

;printMinMax(SIC_month,0)

  SIC(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:,:)=SIC(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:,:)+SIC_month(0:month(imonth)-1,:,:)

delete(SIC_month)

end do
end do

;print(SIC)

S=new((/Nlat0,Nlon/),float)
S=0.
do ilat=0,Nlat0-1
 do ilon=0, Nlon-1
      if (.not.ismissing(SIC(0,ilat,ilon)))
        S(ilat,ilon)=1
end if
end do
end do 

Area=0.
 do ilat=4,9
  do ilon=76,108
 Area=Area+(R*R*cos(lat0(ilat)*0.0174533)*(2.5*0.0174533)*(2.5*0.0174533))*S(ilat,ilon)
end do
end do

;-----------------------------------------------------------
;ENERGY LATENT, PLANETARY WAVES, DAILY RESOLUTION


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

   Tf1=addfile(datapath_energy+"vQtot."+iyx+"."+imx+".WN20.nc","r")

   vQtot=Tf1->vQtot

   vQtot_days=new((/Nwave,month(imonth),Nlat01/),float)
   vQtot_days!0="WaveNumb"
   vQtot_days!1="timeday"
   vQtot_days!2="lat"
   vQtot_days&WaveNumb=Tf1&WaveNumb
   vQtot_days&timeday=fspan(0,month(imonth)-1,month(imonth))
   vQtot_days&lat=Tf1&lat

     do i=0,month(imonth)-1
       vQtot2=vQtot(:,4*i:4*i+3,:)
       vQtot3=dim_sum_n(vQtot2,1)/4
       vQtot_days(:,i,:)=(/vQtot3(:,:)/)
     end do

   vQ_planetary=vQtot_days(1:5,:,:) 
   vQp=dim_sum_n(vQ_planetary,0)
   vQp!0="time"
   vQp!1="lat"

   vQ(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:)=vQ(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:)+vQp(0:month(imonth)-1,:)

delete(vQp)
delete(vQtot)
delete(vQtot_days)
delete(vQtot2)
delete(vQtot3)
delete(vQ_planetary)

end do
end do

;printVarSummary(vQ)

vQp_lon=new((/365*Nyear,Nlat01,6/),float)
vQp_lon!0="time"
vQp_lon!1="lat"
vQp_lon!2="lon"
vQp_lon&time=vQ&time
vQp_lon&lat=Tf1&lat
vQp_lon&lon=fspan(0,5,6)

do ilon=0,5
vQp_lon(:,:,ilon)=vQ(:,:)
end do

do ilat=0,Nlat01-1
 vQp_lon(:,ilat,:)=vQp_lon(:,ilat,:)*2*PI*R*cos(lat(ilat)*2*PI/360.)/1e15
end do

;printVarSummary(vQp_lon)

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

vQp_clim=clmDayTLL(vQp_lon,yyyyddd)         
SIC_clim=clmDayTLL(SIC,yyyyddd)

vQp_sm=runave_n(vQp_clim,8,-1,0)
SIC_sm=runave_n(SIC_clim,8,-1,0)

vQp_an=calcDayAnomTLL(vQp_lon,yyyyddd,vQp_sm)
SIC_an=calcDayAnomTLL(SIC,yyyyddd,SIC_sm)

;printVarSummary(vQp_an)
;printMinMax(vQp_an,0)
;printMinMax(SIC_an,0)  

;-----------------------------------------------------------
;REGRESSION FOR THE WINTER

lag=40
Nlag=2*lag+1


Reg_winter=new((/Nlag,Nlat0,Nlon/),float)
Reg_winter!0="time_lag"
Reg_winter!1="lat"
Reg_winter!2="lon"
Reg_winter&time_lag=fspan(-lag,lag,2*lag+1)
Reg_winter&lat=Tf0&g0_lat_1
Reg_winter&lon=Tf0&g0_lon_2
Reg_winter=0.

do iyear=1,Nyear-2

vQp_reg_winter=vQp_an(343+365*iyear:109+365+365*iyear,40,0) ;70°N and no longitude and summer
SIC_reg_winter=SIC_an(343+365*iyear:109+365+365*iyear,:,:)


  Reg=new((/Nlag,Nlat0,Nlon/),float)
  Reg!0="time_lag"
  Reg!1="lat"
  Reg!2="lon"
  Reg&time_lag=fspan(-lag,lag,2*lag+1)
  Reg&lat=Tf0&g0_lat_1
  Reg&lon=Tf0&g0_lon_2
 ; Reg=0.

  N=dimsizes(vQp_reg_winter&time)

  do ilat=0,Nlat0-1
   do ilon=0,Nlon-1
    do ilag=-lag,lag

  Reg(lag+ilag,ilat,ilon)=avg(SIC_reg_winter(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg_winter(lag:N-1-lag))/(variance(vQp_reg_winter(lag:N-1-lag)))


    end do
   end do
  end do

  Reg_winter=Reg_winter+Reg

end do


;printMinMax(Reg,0)
;printVarSummary(Reg_summer)
;print(Reg&lat)
;print(Reg&lon)

Reg_winter=Reg_winter/(Nyear-2)

printMinMax(Reg_winter,0)


delete(N)


;---------------------------------------------------------
;BARENTZ AREA AND MEAN

Reg_box=Reg_winter(:,4:9,76:108)
Reg_plot0=dim_avg(Reg_box)

lat_barentz=(/80.,77.5,75.,72.5,70.,67.5/)
 do ilat=0,5
Reg_plot0(:,ilat)=Reg_plot0(:,ilat)*cos(lat_barentz(ilat)*0.0174533)
 end do
Reg_plot=dim_sum(Reg_plot0)/Sum_cos_lat

Reg_plot!0="time_lag"
Reg_plot&time_lag=fspan(-lag,lag,2*lag+1)

printMinMax(Reg_plot,0)


;------------------------------------------------------------
;STATISTICAL SIGNIFICANCE FOR WINTER

;RANDOM PHASE FOR vQp_reg


Stat=new((/Nlag/),float)
Stat!0="time_lag"
Stat&time_lag=fspan(-lag,lag,2*lag+1)
Stat=0.

vQp_reg_rand_phase=new(Nyear*365,float)
vQp_reg_rand_phase!0="time"
vQp_reg_rand_phase&time=fspan(0,365*Nyear-1,Nyear*365)


;file_back = "../REGRESSION/STAT_SIGN/Stat_sign_sea_ice.nc"
;a=addfile("../REGRESSION/STAT_SIGN/Stat_sign.nc","r")
;a = addfile(file_back,"r")
;Stat_sign = a->$"Stat_sign"$ ; time,lat,lon
;iSS= a->$"iSS"$
iSS=0
;print(iSS)


 do istat=0,Nstat-1

Reg_winter_rand=new((/Nlag,Nlat0,Nlon/),float)
Reg_winter_rand!0="time_lag"
Reg_winter_rand!1="lat"
Reg_winter_rand!2="lon"
Reg_winter_rand&time_lag=fspan(-lag,lag,2*lag+1)
Reg_winter_rand&lat=Tf0&g0_lat_1
Reg_winter_rand&lon=Tf0&g0_lon_2
Reg_winter_rand=0.


Stat_sign_winter=new((/Nlag/),float)
copy_VarCoords(Stat,Stat_sign_winter)
Stat_sign_winter=0. 

 x=generate_unique_indices(365*Nyear)
 do k=0,365*Nyear-1
  vQp_reg_rand_phase(k)=vQp_an(x(k),40,0)
 end do

do iyear=1,Nyear-2

vQp_reg_winter2=vQp_reg_rand_phase(343+365*iyear:109+365+365*iyear)
SIC_reg_winter2=SIC_an(343+365*iyear:109+365+365*iyear,:,:)

;RGRESSION WITH vQp_reg_rand_phase

  Reg_rand_phase=new((/Nlag,Nlat0,Nlon/),float)
  copy_VarCoords(Reg,Reg_rand_phase)


  N=dimsizes(vQp_reg_winter2&time)

  do ilat=0,Nlat0-1
   do ilon=0,Nlon-1
    do ilag=-lag,lag
  Reg_rand_phase(lag+ilag,ilat,ilon)=avg(SIC_reg_winter2(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg_winter2(lag:N-1-lag))/(variance(vQp_reg_winter2(lag:N-1-lag)))

    end do
   end do
  end do

;printVarSummary(Reg_rand_phase)

Reg_winter_rand=Reg_rand_phase+Reg_winter_rand
delete(Reg_rand_phase)
end do

Reg_winter_rand=Reg_winter_rand/(Nyear-2)


 Reg_rand_phase2=dim_avg(Reg_winter_rand(:,4:9,76:108))    
; Reg_rand_phase_plot=dim_avg(Reg_rand_phase2) 
 do ilat=0,5
Reg_rand_phase2(:,ilat)=Reg_rand_phase2(:,ilat)*cos(lat_barentz(ilat)*0.0174533)
 end do
Reg_rand_phase_plot=dim_sum(Reg_rand_phase2)/Sum_cos_lat


;COMPARISON

  Stat_sign_winter=Stat_sign_winter+where(abs(Reg_rand_phase_plot).lt.abs(Reg_plot),1.,0.)

  Stat=Stat+Stat_sign_winter

print(istat)

delete(Stat_sign_winter)
delete(Reg_winter_rand)

end do

Stat=Stat/Nstat*100

printMinMax(Stat,0)


;  system("rm -f "+file_back)
;  ncdfRaw = addfile(file_back ,"c")
;  fAtt = True
;  fAtt@title = "MC significance test"
;  fAtt@Conventions = "None"
;  fAtt@creation_date = systemfunc("date")
;  fileattdef( ncdfRaw, fAtt )
;  ncdfRaw->Stat_sign = Stat_sign
;  ncdfRaw->iSS = i+1


;---------------------------------------------------------
;PLOT

Stat_sign_plot=Stat
Stat_sign_plot!0="time_lag"
Stat_sign_plot&time_lag=fspan(-lag,lag,2*lag+1)

printMinMax(Reg_plot,0)
printMinMax(Stat_sign_plot,0)


Reg_split=new((/3,Nlag/),float)

Reg_split(0,:)=Reg_plot*Area
Reg_split(1,:)=where(Stat_sign_plot .gt. 95,Reg_plot*Area,Reg_plot@_FillValue)
;Reg_split(2,:)=where(Stat_sign_plot .gt. 99,Reg_plot*Area,Reg_plot@_FillValue)
Reg_split(2,:)=0.

;----------------------------------------------------------
;PLOT 2D


wks  = gsn_open_wks("eps","regression_sea_ice_barentz_winter_PW")
 

res = True 
res@tiMainString="Winter/Latent/PW" 
res@tiXAxisString = "Time lag (days)"  
res@tiYAxisString   = "Regression coefficient (/PW)"  
res@xyLineColors=(/"black","green","black"/)
res@xyLineThicknesses=(/1,4,1/)
res@xyDashPatterns=(/0,0,1/)
     
plot=gsn_csm_xy(wks,Reg_plot&time_lag, Reg_split,res)   

;-------------------------------------------------------------
;SAVE

path_saving = "../REGRESSION/STAT_SIGN/"
system("/bin/rm -f "+path_saving+"Sea_ice_barents_winter_PW.nc")
nc = addfile(path_saving + "Sea_ice_barents_winter_PW.nc","c")
nc->Reg_split_barents_winter_PW=Reg_split
;nc->Stat_winter=Stat_winter





;-----------------------------------------------------------------
print("REGRESSION SEA ICE DONE BARENTZ WINTER PW")

end



