
begin

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;-----------------------------------------------------------
;DEFINE VARIABLE

datapath_temperature="/global/work/runegg/interim/surface/2.5x2.5/an/SAT/"
datapath_energy="/global/work/runegg/interim/EnergySplit/Waves/"
Tf0=addfile(datapath_temperature+"SAT.1980.01.grb","r")
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

T=new((/365*Nyear,Nlat0,Nlon/),float)
T!0="time_year"
T!1="lat"
T!2="lon"
T&time_year=fspan(0,364*Nyear,365*Nyear)
T&lat=Tf0&g0_lat_1
T&lon=Tf0&g0_lon_2
T=0.

vQ=new((/365*Nyear,Nlat01/),float)
vQ!0="time_year"
vQ!1="lat"
vQ&time_year=fspan(0,364*Nyear,365*Nyear)
vQ&lat=Tf01&lat
vQ=0.

delete(Tf01)

;-----------------------------------------------------------
;READ TEMPERATURE


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


  datapath_temperature="/global/work/runegg/interim/surface/2.5x2.5/an/SAT/"
  datapath2_temperature="../nc/"

  fileogrb = datapath_temperature+"SAT."+iyx+"."+imx+".grb"
  fileonc = datapath2_temperature+"SAT."+iyx+"."+imx+".nc"
  system("cdo -R -f nc -t ecmwf copy "+fileogrb+" "+fileonc)
  system("/bin/rm -f "+fileonc)

  Tf2=addfile(datapath_temperature+"SAT."+iyx+"."+imx+".grb","r")
  T_month=Tf2->2T_GDS0_SFC

;print(dimsizes(Tf2&forecast_time0))

  T(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:,:)=T(monthb(imonth)+iyear*365:monthe(imonth)+iyear*365,:,:)+T_month(0:month(imonth)-1,:,:)

delete(T_month)

end do
end do

;printVarSummary(T)


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
vQp_lon&time=vQ&time_year
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
SAT_clim=clmDayTLL(T,yyyyddd)

vQp_sm=runave_n(vQp_clim,8,-1,0)
SAT_sm=runave_n(SAT_clim,8,-1,0)

vQp_an=calcDayAnomTLL(vQp_lon,yyyyddd,vQp_sm)
SAT_an=calcDayAnomTLL(T,yyyyddd,SAT_sm)

;vQp_sm=runave_n(vQp_an,7,-1,0)
;SAT_sm=runave_n(SAT_an,7,-1,0)

;vQp_sm!0="time"
;vQp_sm&time=fspan(0,365*Nyear-1,365*Nyear)
;SAT_sm!0="time"
;SAT_sm&time=fspan(0,Nyear*365-1,Nyear*365)

;printVarSummary(vQp_sm)
;printVarSummary(vQp_an)
;printMinMax(vQp_an,0)
;printMinMax(SAT_an,0)

;-----------------------------------------------------------
;REGRESSION_SUMMER

lag=20
Nlag=2*lag+1

Reg_summer=new((/Nlag,Nlat0,Nlon/),float)
Reg_summer!0="time_lag"
Reg_summer!1="lat"
Reg_summer!2="lon"
Reg_summer&time_lag=fspan(-lag,lag,2*lag+1)
Reg_summer&lat=Tf0&g0_lat_1
Reg_summer&lon=Tf0&g0_lon_2
Reg_summer=0.
Reg_winter=new((/Nlag,Nlat0,Nlon/),float)
copy_VarCoords(Reg_summer,Reg_winter)
Reg_winter=0.

do iyear=0,Nyear-1

  vQp_reg_summer=vQp_an(161+365*iyear:275+365*iyear,40,0) ;70°N and no longitude and summer
  SAT_reg_summer=SAT_an(161+365*iyear:275+365*iyear,:,:)


  Reg=new((/Nlag,Nlat0,Nlon/),float)
  Reg!0="time_lag"
  Reg!1="lat"
  Reg!2="lon"
  Reg&time_lag=fspan(-lag,lag,2*lag+1)
  Reg&lat=Tf0&g0_lat_1
  Reg&lon=Tf0&g0_lon_2

  N=dimsizes(vQp_reg_summer&time_year)

   do ilat=0,Nlat0-1
    do ilon=0,Nlon-1
     do ilag=-lag,lag
  Reg(lag+ilag,ilat,ilon)=avg(SAT_reg_summer(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg_summer(lag:N-1-lag))/(variance(vQp_reg_summer(lag:N-1-lag)))

  end do
  end do
  end do

  Reg_summer=Reg_summer+Reg  

end do

;-------------------------------------------------------------
;REGRESSION_WINTER

do iyear=1,Nyear-2

  vQp_reg_winter=vQp_an(343+365*iyear:109+365+365*iyear,40,0) ;70°N and no longitude and summer
  SAT_reg_winter=SAT_an(343+365*iyear:109+365+365*iyear,:,:)
;printVarSummary(vQp_reg)
;printVarSummary(SAT_reg)

  Reg2=new((/Nlag,Nlat0,Nlon/),float)
  Reg2!0="time_lag"
  Reg2!1="lat"
  Reg2!2="lon"
  Reg2&time_lag=fspan(-lag,lag,2*lag+1)
  Reg2&lat=Tf0&g0_lat_1
  Reg2&lon=Tf0&g0_lon_2

  N=dimsizes(vQp_reg_winter&time_year)

   do ilat=0,Nlat0-1
    do ilon=0,Nlon-1
     do ilag=-lag,lag
  Reg2(lag+ilag,ilat,ilon)=avg(SAT_reg_winter(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg_winter(lag:N-1-lag))/(variance(vQp_reg_winter(lag:N-1-lag)))

  end do
  end do
  end do

  Reg_winter=Reg_winter+Reg2  

end do

Reg_summer=Reg_summer/Nyear
Reg_winter=Reg_winter/(Nyear-2)


;printMinMax(Reg_summer,0)
;printVarSummary(Reg_summer)


;------------------------------------------------------------
;STATISTICAL SIGNIFICANCE FOR SUMMER

;RANDOM PHASE FOR vQ_reg


vQp_reg_rand_phase_summer=new(115,float)
vQp_reg_rand_phase_summer!0="time"
vQp_reg_rand_phase_summer&time=fspan(0,114,115)

Stat_sign_summer=new((/Nlag,Nlat0/),float)
Stat_sign_summer!0="time_lag"
Stat_sign_summer!1="lat"
Stat_sign_summer&time_lag=fspan(-lag,lag,2*lag+1)
Stat_sign_summer&lat=Tf0&g0_lat_1
Stat_sign_summer=0.

Stat_summer=new((/Nlag,Nlat0/),float)
copy_VarCoords(Stat_sign_summer,Stat_summer)
Stat_summer=0.

;file_back = "../REGRESSION/STAT_SIGN/Stat_sign_summer.nc"
;a=addfile("../REGRESSION/STAT_SIGN/Stat_sign_summer.nc","r")
;a = addfile(file_back,"r")
;Stat_sign = a->$"Stat_summer"$ ; time,lat,lon
;iSS= a->$"isummer"$
isummer=-1



;printMinMax(Stat_sign,0)

;do iyear=isummer+1,Nyear-1

Stat_sign_summer=new((/Nlag,Nlat0/),float)
Stat_sign_summer!0="time_lag"
Stat_sign_summer!1="lat"
Stat_sign_summer&time_lag=fspan(-lag,lag,2*lag+1)
Stat_sign_summer&lat=Tf0&g0_lat_1
Stat_sign_summer=0.
;vQp_reg_summer2=vQp_an(161+365*iyear:275+365*iyear,40,0)
;SAT_reg_summer2=SAT_an(161+365*iyear:275+365*iyear,:,:)

;print(vQp_reg_summer)

;  do i=0,Nstat-1

; x=generate_unique_indices(115)
; do k=0,114
;  vQp_reg_rand_phase_summer(k)=vQp_reg_summer2(x(k))
; end do


;RGRESSION WITH vQp_reg_rand_phase_summer

;  Reg_rand_phase_summer=new((/Nlag,Nlat0,Nlon/),float)
;  copy_VarCoords(Reg_summer,Reg_rand_phase_summer)


;  N=dimsizes(vQp_reg_rand_phase_summer&time_year)

;  do ilat=0,Nlat0-1
;   do ilon=0,Nlon-1
;    do ilag=-lag,lag
;  Reg_rand_phase_summer(lag+ilag,ilat,ilon)=avg(SAT_reg_summer2(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg_rand_phase_summer(lag:N-1-lag))/(variance(vQp_reg_rand_phase_summer(lag:N-1-lag)))

;    end do
;   end do
;  end do

;printVarSummary(Reg_rand_phase_summer)



;COMPARISON

 ; Stat_sign_summer=Stat_sign_summer+where(abs(dim_avg(Reg_rand_phase_summer)).lt.abs(dim_avg(Reg_summer)),1.,0.)


;end do

;Stat_sign_summer=Stat_sign_summer/Nstat*100

;printMinMax(Stat_sign_summer,0)


;Stat_summer=Stat_summer+Stat_sign_summer

;printMinMax(Stat_summer,0)
;delete(Stat_sign_summer)
;print(iyear)

 ;system("rm -f "+file_back)
 ;ncdfRaw = addfile(file_back ,"c")
 ;fAtt = True
 ;fAtt@title = "MC significance test"
 ;fAtt@Conventions = "None"
 ;fAtt@creation_date = systemfunc("date")
 ;fileattdef( ncdfRaw, fAtt )

 ;ncdfRaw->Stat_summer = Stat_summer
 ;ncdfRaw->isummer = i


;end do

;Stat_summer=Stat_summer/Nyear


;printMinMax(Stat_summer,0)


;------------------------------------------------------------
;STATISTICAL SIGNIFICANCE FOR WINTER

;RANDOM PHASE FOR vQ_reg

Stat=new((/Nlag,25/),float)
Stat!0="time_lag"
Stat!1="lat"
Stat&time_lag=fspan(-lag,lag,2*lag+1)
Stat&lat=fspan(90,30,25)
Stat=0.

;file_back = "../REGRESSION/STAT_SIGN/Stat_sign_winter.nc"
;a=addfile("../REGRESSION/STAT_SIGN/Stat_sign_winter.nc","r")
;a = addfile(file_back,"r")
;Stat_sign = a->$"Stat_winter"$ ; time,lat,lon
;iwinter= a->$"iwinter"$
iwinter=1


vQp_reg_rand_phase_winter=new(365*Nyear,float)
vQp_reg_rand_phase_winter!0="time_year"
vQp_reg_rand_phase_winter&time_year=fspan(0,365*Nyear-1,365*Nyear)
;copy_VarCoords(vQp_reg_rand_phase_winter,vQp_reg_rand_phase_winter)

do istat=0,Nstat-1

Stat_winter=new((/Nlag,25/),float)
Stat_winter!0="time_lag"
Stat_winter!1="lat"
Stat_winter&time_lag=fspan(-lag,lag,2*lag+1)
Stat_winter&lat=fspan(90,30,25)
Stat_winter=0.

 y=generate_unique_indices(365*Nyear)
 do k=0,365*Nyear-1
  vQp_reg_rand_phase_winter(k)=vQp_an(y(k),40,0)
 end do

do iyear=iwinter,Nyear-2
 
 Stat_sign_winter=new((/Nlag,25/),float)
 Stat_sign_winter!0="time_lag"
 Stat_sign_winter!1="lat"
 Stat_sign_winter&time_lag=fspan(-lag,lag,2*lag+1)
 Stat_sign_winter&lat=fspan(90,30,25)
 Stat_sign_winter=0.

 vQp_reg_winter2=vQp_reg_rand_phase_winter(343+365*iyear:109+365+365*iyear)
 SAT_reg_winter2=SAT_an(343+365*iyear:109+365+365*iyear,:,:)

; y=generate_unique_indices(132)
; do k=0,131
;  vQp_reg_rand_phase_winter(k)=vQp_reg_winter2(y(k))
; end do

;print(vQp_reg_rand_phase_winter(5))

;RGRESSION WITH vQp_reg_rand_phase_summer

  Reg_rand_phase_winter=new((/Nlag,Nlat0,Nlon/),float)
  copy_VarCoords(Reg_winter,Reg_rand_phase_winter)


  N=dimsizes(vQp_reg_winter2&time_year)

  do ilat=0,Nlat0-1
   do ilon=0,Nlon-1
    do ilag=-lag,lag
  Reg_rand_phase_winter(lag+ilag,ilat,ilon)=avg(SAT_reg_winter2(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg_winter2(lag:N-1-lag))/(variance(vQp_reg_winter2(lag:N-1-lag)))

    end do
   end do
  end do

;printVarSummary(Reg_rand_phase_winter)

;COMPARISON

  Stat_sign_winter=Stat_sign_winter+where(abs(dim_avg(Reg_rand_phase_winter(:,0:24,:))).lt.abs(dim_avg(Reg_winter(:,0:24,:))),1.,0.)


Stat_winter=Stat_winter+Stat_sign_winter
delete(Stat_sign_winter)

end do

Stat_winter=Stat_winter/(Nyear-2)

Stat=Stat+Stat_winter
delete(Stat_winter)
delete(vQp_reg_rand_phase_winter)
print(istat)
end do

Stat=Stat/Nstat*100


printMinMax(Stat,0)






 



; system("rm -f "+file_back)
; ncdfRaw = addfile(file_back ,"c")
; fAtt = True
; fAtt@title = "MC significance test"
; fAtt@Conventions = "None"
; fAtt@creation_date = systemfunc("date")
; fileattdef( ncdfRaw, fAtt )
; ncdfRaw->Stat_winter = Stat_winter
; ncdfRaw->iwinter = i+1





;-----------------------------------------------------------
;PLOT

Reg2_summer=dim_avg(Reg_summer(:,0:24,:))
Reg2_summer!0="time_lag"
Reg2_summer!1="lat"
Reg2_summer&time_lag=fspan(-20,20,41)
Reg2_summer&lat=fspan(90,30,25)

Reg2_winter=dim_avg(Reg_winter(:,0:24,:))
copy_VarCoords(Reg2_summer,Reg2_winter)

printMinMax(Reg2_summer,0)
printMinMax(Reg2_winter,0)

Stat_summer2=Stat_summer(:,0:24)
;Stat_winter2=Stat_winter(:,0:24)

printMinMax(Stat_summer,0)
printMinMax(Stat,0)
;-------------------------------------------------------------
;DRAW

wks  = gsn_open_wks("eps","regression_summer_winter")
plot=new(2,graphic)
res = True         
res@gsnDraw=False
res@gsnFrame=False
res@cnInfoLabelOn=False

res@cnFillOn = True       
res@cnFillPalette = "BlWhRe" 
;res@vpWidthF = 0.2          
;res@vpHeightF = 0.1
res@gsnMaximize = True
res@lbLabelBarOn = False 
res@cnLineLabelsOn=False
res@cnLevelSelectionMode ="ManualLevels"
res@cnMaxLevelValF = 7.
res@cnMinLevelValF = -7.
res@cnLevelSpacingF = 1
res@tiXAxisString = "latitude (deg)"
res@tiYAxisString = "time lag"
res@trXReverse=True

plot1 = gsn_csm_contour(wks,Reg2_summer,res) 

res2                      = True
res2@gsnDraw=False
res2@gsnFrame=False 
res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
res2@cnMinLevelValF       = 95.              ; set min contour level
res2@cnMaxLevelValF       = 99.              ; set max contour level
res2@cnLevelSpacingF      = 4.              ; set contour spacing
res2@cnLineLabelsOn       = False
res2@cnMonoLineColor      = False
res2@cnLineColors          = (/"Green","Green"/)
res2@cnLineThicknessF     = 2.
res2@cnInfoLabelOn        = False
res2@cnMonoLineDashPattern = False
res2@cnLineDashPatterns   = (/2,0/)
res2@tiXAxisString = "latitude (deg)"
res2@tiYAxisString = "time lag"
res2@trXReverse=True

;plot2=gsn_csm_contour(wks,Stat_summer2,res2)
plot3 = gsn_csm_contour(wks,Reg2_winter,res) 
plot4=gsn_csm_contour(wks,Stat,res2)


;overlay(plot1,plot2)
overlay(plot3,plot4)

plot(0)=plot1
plot(1)=plot3

          
resP=True
resP@gsnFrame=False
resP@txFontThicknessF=2
resP@txFontHeightF=0.02
resP@txString="Summer and winter regression (K/PW)"
resP@lbLabelFontHeightF  = 0.007
resP@gsnPanelLabelBar=True
;resP@tiXAxisString = "latitude (deg)"
;resP@tiYAxisString = "time lag"


gsn_panel(wks,plot,(/2,1/),resP)   
frame(wks)
;-----------------------------------------------------------
;SAVE

path_saving = "../REGRESSION/STAT_SIGN/"
system("/bin/rm -f "+path_saving+"Stat_sign_season.nc")
nc = addfile(path_saving + "Stat_sign_season.nc","c")
nc->Stat_winter=Stat
;nc->Stat_summer=Stat_summer





print("REGREESION SEASON DONE")

end

