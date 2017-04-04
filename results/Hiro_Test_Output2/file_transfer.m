prof = [45];
tide = [1 8];
height = [1 3];
resi = [1e-1 1e-0];
wea  = [1 3 5];
wea_v= [1e-4 1e-3 1e-2]
co   = [1 3 5];

in_path = 'S:\env\Share\Hiro\wvsw\20170403\';
ou_path = 'C:\Users\hmat258\Desktop\Martin test\';

for pa=1:length(prof)
for ta=1:length(tide)
for ha=1:length(height)
for ra=1:length(resi)
for wa=1:length(wea)
for ca=1:length(co)
    
    file_path = [in_path,num2str(prof(pa)),'\tide-',num2str(tide(ta)),'\pressure1\wea1\weasd',num2str(wea(wa)),'\co',num2str(co(ca)),'\'];
    file      = ['p',num2str(prof(pa)),'-',num2str(height(ha)),'m-resi',num2str(resi(ra)),'-save_prof.txt'];
    new_name = [ou_path,'P',num2str(prof(pa)),'-T',num2str(tide(ta)),'m-wea',num2str(wea(wa)),'-co',num2str(co(ca)),'-wave',num2str(height(ha)),'m-resi',num2str(resi(ra)*0.1),'-save_prof.txt'];
    
    copyfile([file_path,file],new_name);
end
end
end
end
end
end
    
    
    
  
    
