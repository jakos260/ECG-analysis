bouwen van de windows 64 bits:
gebruik visual studio 10++ (er is een gratis expres versie)
verander eerst
#IFDEF WIN32 
to 
#IFDEF _WIN64
en dan met dit commando
mex pnet.c ws2_32.lib

