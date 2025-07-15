clear all
% Matlab client for calling ECE module of TRAVIS-code 
%----------------------------------------------------------

%----------------------------------------------------------
% --  travisECESpectrum example

% -- set path to java-jar file.
javaaddpath H:\W7XTravis\matlab\Travis.jar

import ipp.w7x.travis.service.*;
serv=Travis;
%serv=Travis(java.net.URL('http://127.0.0.1:55221/Travis?wsdl'));

port=serv.getTravis;
methodsview TravisPortType; 

input = ECEinputTravis;

me = MagneticEquilibriumTravis;
ac = AntennaConfigurationTravis;
pr = PlasmaProfilesTravis;
fm = WaveFrequencyModeTravis;


me.setUrl('http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/run/w7x.1000_1000_1000_1000_+0500_+0500.06.088/wout.txt');
me.setB0Type('at angle on magn.axis');
me.setPhiRef( java.lang.Double(0));
me.setB0Ref( java.lang.Double(2.5));

ac.setAntennaCoordType('cyl'); 
position = ac.getAntennaPosition();
position.add(java.lang.Double(6.5));
position.add(java.lang.Double(6.1));
position.add(java.lang.Double(0.35));

ac.setTargetPositionType('cyl'); 
tPosition = ac.getTargetPosition();
tPosition.add(java.lang.Double(5.5480));
tPosition.add(java.lang.Double(5.2392));
tPosition.add(java.lang.Double(0.0576));

ac.setRbeam(java.lang.Double(0.02));
ac.setRfocus(java.lang.Double(1.0));
ac.setQOemul(java.lang.Integer(0));


pr.setLabelType('tor_rho');
rho = pr.getRho();
ne = pr.getNe();
te = pr.getTe();

n = 50;
for i = 0:n 
    x = i /double(n-1);
    rho.add(java.lang.Double(x));
    ne.add(java.lang.Double(2e19*(1-power(x, 2))));
    te.add(java.lang.Double(5*(1-power(x, 2))));
end


nf = 21;
df = 4;
frequencies = fm.getFrequencies();
waveModeXorO = fm.getWaveMode();

for i = 0:nf
  frequencies.add(java.lang.Double(100+i*df));
  waveModeXorO.add(java.lang.Integer(1));
end

input.setMagneticEquilibrium(me);
input.setAntennaConfiguration(ac);
input.setPlasmaProfiles(pr);
input.setWaveFreqMode(fm);

eceinp = TravisECESpectrum;
fact = ObjectFactory; 
inp1 = fact.createTravisECESpectrumInput(input);
eceinp.setInput(inp1);


% Start travis call....
response = port.travisECESpectrum(eceinp); 
% Ready

Tece = response.getTece();

f=frequencies.iterator;
t=Tece.iterator;

figure;
hold on;
while f.hasNext > 0
    plot(f.next, t.next, '-*b');
end
hold off
%pause;


%==================================================================
% --  travisECESpectrumAll example

eceinp2 = TravisECESpectrumAll;
fact = ObjectFactory; 
inp2 = fact.createTravisECESpectrumAllInput(input);
eceinp2.setInput(inp2);


% Start travis call....
response = port.travisECESpectrumAll(eceinp2); 
% Ready

Tece = response.getTece();
Rcy = response.getRcy();


figure;
hold on;
f=frequencies.iterator;
t=Tece.iterator;
while f.hasNext > 0
    plot(f.next, t.next, '-*b');
end
hold off
%pause;

rho = response.getRhoeceb();
figure;
hold on;
r=rho.iterator;
t=Tece.iterator;
while t.hasNext > 0
    rl=r.next;   % left
    r0=r.next;   % central
    rr=r.next;   % right
    plot(r0, t.next, '-*b');
end
hold off
%pause;


%list of cached configurations
l=ListCache;
l.setPattern('s')
res=port.listCache(l);
res.getLst();
ans;

%==================================================================
% GetEpsEff example         
%methodsview GetEpsEff;   
%methodsview Eeff;

inp = GetEpsEff;
inp.setUrl('http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/run/w7x.1000_1000_1000_1000_+0500_+0500.06.088/wout.txt');
resp=port.getEpsEff(inp);
e=resp.getE();
s=resp.getS();


figure;
hold on;
i=s.iterator;
j=e.iterator;
while i.hasNext > 0
    plot(sqrt(i.next), j.next, '-*r','LineWidth',2);
end
hold off
%pause;


