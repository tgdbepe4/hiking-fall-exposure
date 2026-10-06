# Erzeugt synthetisches Testgelaende (dem.tif, forest_4326.tif, ways.gpkg) im aktuellen Ordner.
# Erwartete Klassen: A_flach=1, B_35Grad=2, C_45Grad=3, E_45Grad_Wald=3, F_Felskante=4, G_Plateau=1
import numpy as np, rasterio, geopandas as gpd
from rasterio.transform import from_origin
from rasterio.warp import reproject, calculate_default_transform, Resampling
from shapely.geometry import LineString
X0,Y0,res=2600000,1200000,0.5
W,H=4000,2000
xs=np.arange(W)*res+res/2; ys=1000-(np.arange(H)*res+res/2)
X,Y=np.meshgrid(xs,ys)
z=np.full(X.shape,500.0)
t=lambda d: np.tan(np.radians(d))
r=X<500; z[r]=500+Y[r]*t(5)
r=(X>=500)&(X<900); z[r]=400+Y[r]*t(35)
r=(X>=900)&(X<1300); z[r]=300+Y[r]*t(45)
r=X>=1300
zc=np.where(Y>=500,700.0,np.where(Y>=497,650+(Y-497)/3*50,650-(497-Y)*t(20)))
z[r]=zc[r]
z+=np.random.default_rng(1).normal(0,0.05,z.shape)
prof=dict(driver="GTiff",width=W,height=H,count=1,dtype="float32",crs="EPSG:2056",
          transform=from_origin(X0,Y0+1000,res,res),nodata=-9999)
with rasterio.open("dem.tif","w",**prof) as d: d.write(z.astype("float32"),1)
# Wald 10 m, 100 % fuer x 1100..1300, sonst 0 -> in EPSG:4326 speichern (Reprojektionstest)
fw,fh=200,100
fx=np.arange(fw)*10+5
f=np.zeros((fh,fw),"float32"); f[:,(fx>=1100)&(fx<1300)]=100
st=from_origin(X0,Y0+1000,10,10)
dt,dw,dh=calculate_default_transform("EPSG:2056","EPSG:4326",fw,fh,X0,Y0,X0+2000,Y0+1000)
dst=np.full((dh,dw),255,"float32")
reproject(f,dst,src_transform=st,src_crs="EPSG:2056",dst_transform=dt,dst_crs="EPSG:4326",
          resampling=Resampling.nearest,dst_nodata=255)
with rasterio.open("forest_4326.tif","w",driver="GTiff",width=dw,height=dh,count=1,
     dtype="float32",crs="EPSG:4326",transform=dt,nodata=255) as d: d.write(dst,1)
L=lambda x0,x1,y: LineString([(X0+x0,Y0+y),(X0+x1,Y0+y)])
names=["A_flach","B_35Grad","C_45Grad","E_45Grad_Wald","F_Felskante","G_Plateau"]
geoms=[L(100,450,500),L(550,850,500),L(950,1080,500),L(1130,1270,500),L(1400,1900,503),L(1400,1900,700)]
gpd.GeoDataFrame({"osm_id":range(1,7),"name":names,"sac_scale":"mountain_hiking"},
  geometry=geoms,crs="EPSG:2056").to_crs(4326).to_file("ways.gpkg",layer="sac_scale_test")
