function [c,P] = myTransformada(c,P)
    fs=50000;
    ts=1/fs;
    t=-1/50:ts:1/50;
    
    M=0;
    n=0;
    w=(-1000*P:1000*P)*2*pi;
    m=c;
    
    for tt=t
        n=n+1;
        M=M+m(n)*exp(-j*w*tt)*ts;
    end
    k=w/(2*pi);
    d=abs (M);
    plot (k,d,'k','LineWidth',1.3);
    
end

