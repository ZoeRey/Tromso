load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

path="../ARCTICLE/BACK-UP/"

Tf1 = addfile(path+"vQ.nc","r")
Tf2 = addfile(path+"vE.nc","r")
Tf3 = addfile(path+"Etot.nc","r")

Etot= Tf3->Etot
vQ = Tf1->vQ
vE = Tf2->vE

data=new((/3,dimsizes(Etot&lat)/),float)

data(0,:)=Etot
data(1,:)=vQ
data(2,:)=vE

wks = gsn_open_wks( "eps", "Total+DS+L")
res=True
res@tiMainString="Energy Transport in PW as a function of lat"

res@xyLineThicknesses = (/  5.0,   1.0,   1.0/)
res@xyLineColors      = (/"red","blue","orange"/)


plot = gsn_csm_xy (wks, Tf3&lat, data, res)

end




