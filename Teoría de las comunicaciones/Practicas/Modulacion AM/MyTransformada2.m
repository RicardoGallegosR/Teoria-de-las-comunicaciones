function [ fx,N] = MyTransformada2( fx,N)
    fs=2500;
    ts=1/fs;
    m=fx;
    M= fftshift(fft(m,N))*ts;
    f=linspace(-fs/2, fs/2, N);
    plot(f, abs(M),'k','LineWidth',1.1);
    grid on
end
