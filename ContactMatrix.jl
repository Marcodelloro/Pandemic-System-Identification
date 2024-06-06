# For info on age groups https://covid19.infn.it/iss/
import XLSX
using LinearAlgebra

c_ovr = XLSX.readdata("ContactMarticesRaw.xlsx", "Overall!A1:P16")
c_home = XLSX.readdata("ContactMarticesRaw.xlsx", "Home!A1:P16")
c_sch= XLSX.readdata("ContactMarticesRaw.xlsx", "School!A1:P16")
c_work = XLSX.readdata("ContactMarticesRaw.xlsx", "Work!A1:P16")
c_othr = XLSX.readdata("ContactMarticesRaw.xlsx", "Others!A1:P16")

x = [ones(1,8) zeros(1,8); 
    zeros(1,8) ones(1,4) zeros(1,4);
    zeros(1,8) zeros(1,4) ones(1,4)]

# 3x3 matrix rearranging to match available data
c_ovr = x*c_ovr*x' 
c_home = x*c_home*x'
c_sch= x*c_sch*x'
c_work = x*c_work*x'
c_othr =x*c_othr*x'