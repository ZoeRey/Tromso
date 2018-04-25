


begin

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


;a=addfile("../REGRESSION/STAT_SIGN/Stat_sign.nc","r")
;Stat_sign = a->$"Stat_sign"$ ; time,lat,lon
;printMinMax(Stat_sign,0)

;file_back = "../REGRESSION/STAT_SIGN/6.nc"
;a=addfile("../REGRESSION/STAT_SIGN/Stat_sign.nc","r")
;a = addfile(file_back,"r")
;Stat_sign = a->$"Stat_sign"$ ; time,lat,lon
;printMinMax(Stat_sign,0)

;-----------------------------------------------------------
;DEFINE VARIABLE

datapath_temperature="/global/work/runegg/interim/surface/2.5x2.5/an/SAT/"
datapath_energy="/global/work/runegg/interim/EnergySplit/Waves/"
Tf0=addfile(datapath_temperature+"SAT.1980.01.grb","r")
Tf01=addfile(datapath_energy+"vQtot.1980.01.WN20.nc","r")

PI=3.14159
R=6371000.

;THE NUMBER OF YEARS MUST BE EVEN IN ORDER TO HAVE AN EVEN NUMBER OF DAYS
Nyear=33
Nmonths=12
Nstat=2000

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
  fileonc = datapath2_temperature+"SAT.."+iyx+"."+imx+".nc"
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

;printVarSummary(vQp_an)
;printMinMax(vQp_an,0)
;printMinMax(SAT_an,0)

;-----------------------------------------------------------
;REGRESSION

vQp_reg=vQp_an(:,40,0) ;70°N and no longitude
;printVarSummary(vQp_reg)
;printVarSummary(SAT_an)

lag=20
Nlag=2*lag+1

Reg=new((/Nlag,Nlat0,Nlon/),float)
Reg!0="time_lag"
Reg!1="lat"
Reg!2="lon"
Reg&time_lag=fspan(-lag,lag,2*lag+1)
Reg&lat=Tf0&g0_lat_1
Reg&lon=Tf0&g0_lon_2

N=dimsizes(vQp_an&time_year)

 do ilat=0,Nlat0-1
  do ilon=0,Nlon-1
   do ilag=-lag,lag
Reg(lag+ilag,ilat,ilon)=avg(SAT_an(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg(lag:N-1-lag))/(variance(vQp_reg(lag:N-1-lag)))

end do
end do
end do

;printMinMax(Reg,0)
;printVarSummary(Reg)
;print(Reg&lat)

Reg_negative=new((/Nlag,Nlat0/),float)
Reg_positive=new((/Nlag,Nlat0/),float)
copy_VarCoords(Reg,Reg_positive)
copy_VarCoords(Reg,Reg_negative)
Reg_negative=0.
Reg_positive=0.

Reg4=dim_avg(Reg)


do ilag=0,2*lag
 do ilat=0,Nlat0-1
  do ilon=0,Nlon-1

    if (Reg4(ilag,ilat) .lt. 0) then
      Reg_negative(ilag,ilat)=Reg4(ilag,ilat)
    else
      Reg_positive(ilag,ilat)=(Reg4(ilag,ilat))
    end if

end do
end do
end do

;------------------------------------------------------------
;STATISTICAL SIGNIFICANCE  

;RANDOM PHASE FOR vQ_reg


Stat_sign_positive=new((/Nlag,Nlat0/),float)
Stat_sign_positive!0="time_lag"
Stat_sign_positive!1="lat"
Stat_sign_positive&time_lag=fspan(-lag,lag,2*lag+1)
Stat_sign_positive&lat=Tf0&g0_lat_1
Stat_sign_positive=0.
Stat_sign_negative=new((/Nlag,Nlat0/),float)
copy_VarCoords(Stat_sign_positive,Stat_sign_negative)
Stat_sign_negative=0.


file_back = "../REGRESSION/STAT_SIGN/Stat_sign.nc"
a=addfile("../REGRESSION/STAT_SIGN/Stat_sign.nc","r")
a = addfile(file_back,"r")
Stat_sign_positive = a->$"Stat_sign_positive"$ ; time,lat,lon
Stat_sign_negative = a->$"Stat_sign_negative"$
iSS= a->$"iSS"$

;iSS=0

;printMinMax(Stat_sign,0)

vQp_reg_rand_phase=new(Nyear*365,float)
vQp_reg_rand_phase!0="time"
vQp_reg_rand_phase&time=fspan(0,Nyear*365-1,Nyear*365)


 do i=iSS,Nstat-1

 x=generate_unique_indices(Nyear*365)
 do k=0,Nyear*365-1
  vQp_reg_rand_phase(k)=vQp_reg(x(k))
 end do


;RGRESSION WITH vQp_reg_rand_phase

  Reg_rand_phase=new((/Nlag,Nlat0,Nlon/),float)
  copy_VarCoords(Reg,Reg_rand_phase)


  N=dimsizes(vQp_reg_rand_phase&time_year)

  do ilat=0,Nlat0-1
   do ilon=0,Nlon-1
    do ilag=-lag,lag
  Reg_rand_phase(lag+ilag,ilat,ilon)=avg(SAT_an(lag+ilag:N-lag-1+ilag,ilat,ilon)*vQp_reg_rand_phase(lag:N-1-lag))/(variance(vQp_reg_rand_phase(lag:N-1-lag)))

    end do
   end do
  end do

;printVarSummary(Reg_rand_phase)

;COMPARISON

  Stat_sign_positive=Stat_sign_positive+where(abs(dim_avg(Reg_rand_phase)).lt.abs((Reg_positive)),1.,0.)
  Stat_sign_negative=Stat_sign_negative+where(abs(dim_avg(Reg_rand_phase)).lt.abs((Reg_negative)),1.,0.)




  system("rm -f "+file_back)
  ncdfRaw = addfile(file_back ,"c")
  fAtt = True
  fAtt@title = "MC significance test"
  fAtt@Conventions = "None"
  fAtt@creation_date = systemfunc("date")
  fileattdef( ncdfRaw, fAtt )

  ncdfRaw->Stat_sign_positive=Stat_sign_positive
  ncdfRaw->Stat_sign_negative=Stat_sign_negative
  ncdfRaw->iSS = i+1



print(i)
end do


Stat_sign_positive=Stat_sign_positive/Nstat*100
Stat_sign_negative=Stat_sign_negative/Nstat*100

;-----------------------------------------------------------
;PLOT

Reg2=dim_avg(Reg(:,0:24,:))
Reg2!0="time_lag"
Reg2!1="lat"
Reg2&time_lag=fspan(-lag,lag,2*lag+1)
Reg2&lat=fspan(90,30,25)


printMinMax(Reg2,0)


Stat_sign_positive2=Stat_sign_positive(:,0:24)
Stat_sign_negative2=Stat_sign_negative(:,0:24)
Stat_sign_positive2!0="time_lag"
Stat_sign_positive2&time_lag=fspan(-lag,lag,2*lag+1)
Stat_sign_positive2!1="lat"
Stat_sign_positive2&lat=fspan(90,30,25)
copy_VarCoords(Stat_sign_negative2,Stat_sign_positive2)



;-------------------------------------------------------------
;DRAW

wks  = gsn_open_wks("eps","regresion+SS")

res = True
res@gsnDraw=False
res@gsnFrame=False  
res@cnFillOn = True       
res@cnFillPalette = "BlWhRe" 
res@vpWidthF = 0.5          
res@vpHeightF = 0.3
res@gsnMaximize = True  
res@cnLevelSelectionMode = "ManualLevels"
res@lbLabelBarOn = True       
res@cnMaxLevelValF = 3
res@cnMinLevelValF = -3
res@cnLevelSpacingF = 0.5 
res@trXReverse=True
res@tiXAxisString = "latitude (deg)"
res@tiYAxisString = "time lag"
res@tiMainString    = "Regression (K/PW)"
plot1 = gsn_csm_contour(wks,Reg2,res)

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


plot2=gsn_csm_contour(wks,Stat_sign_positive2,res2)
plot3=gsn_csm_contour(wks,Stat_sign_negative2,res2)

overlay(plot2,plot3)
overlay(plot1,plot2) 
  

draw(plot1)
frame(wks)
;-----------------------------------------------------------
;SAVE
;path_saving = "../REGRESSION/STAT_SIGN/"
;system("/bin/rm -f "+path_saving+"Stat_sign.nc")
;nc = addfile(path_saving + "Stat_sign.nc","c")
;nc->Stat_sign_positive=Stat_sign_positive
;nc->Stat_sign_negative=Stat_sign_negative

print("REGRESSION DONE")

end