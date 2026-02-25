function [S,E,F, D] = ieee754(x)
%IEEE754 Decompose a double and single precision floating point number.
% [S,E,F] = IEEE754(X) returns the sign bit, exponent, and mantissa of an
% IEEE 754 floating point binary array X, expressed as binary digit strings 
% of length 1, 11, and 52, respectively for 64 bitstreams and 1, 8 and 23 
% for 32 bitstreams. 
% Input X is a binary string (1/0). 
%
% [S,E,F, deccie] = IEEE754(X,'bin') returns S, E, and F as floating-point 
% numbers and deccie of the decoded floating point. 
%
% D is equal to (in exact arithmetic and decimal notation)
%
%      (-1)^S * (1 + F/(2^52)) *  2^(E-1023), (64 bits)
%      (-1)^S * (1 + F/(2^23)) *  2^(E-127),  (32 bits)
%
if ~isreal(x) || ~isa(x,'double')
  error('Real, double input required.')
end

if length(x) == 64
    bitstr = reshape(x',[1 64]);  % string of 64 bits in order
    if nargout<2
      S = bitstr;      
    else
      s = bitstr(1);
      e = bitstr(2:12);
      f = bitstr(13:64);
      S = bin2dec(num2str(s));  
      E = bin2dec(num2str(e));  
      F = bin2dec(num2str(f));
      D = (-1)^S*2^(E-1023)*(1+F/2^52);
    end
elseif length(x) == 32
    bitstr = reshape(x',[1 32]);  % string of 32 bits in order
    if nargout<2
      s = bitstr;      
    else
      s = bitstr(1);
      e = bitstr(2:9);
      f = bitstr(10:32);
      S = bin2dec(num2str(s));  
      E = bin2dec(num2str(e));  
      F = bin2dec(num2str(f));
      D = (-1)^S*2^(E-127)*(1+F/2^23);
    end
end
end