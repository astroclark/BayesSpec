
%% Generate data
B1=10000;
alpha=0.03;
gendata(B1,alpha_inj)

return
%% Evalute spectra
analysis='lorentz';

switch analysis
    case 'sine'
        [logprob, w] = sinespec(data,fs);
    case 'lorentz'
%         [samples,acc]=lorentzspec(data,fs,5000,'mh');
        [samples,acc]=lorentzspec(data,fs,5000);

    figure('Position', [100, 50, 1000, 800]);

    % --- Frequencies
    subplot(3,3,1)
    xi=0:0.0005:0.5;
    f=ksdensity(samples(:,1),xi,'support','positive');
    % f=ksdensity(samples,xi);
    plot(xi,f)
    hold on
    plot([0.3,0.3],get(gca,'Ylim'),'r')
    xlabel('w')

    subplot(3,3,2)
    hist(samples(:,1),500)
    hold on
    plot([0.3,0.3],get(gca,'Ylim'),'r')
    xlabel('w')

    subplot(3,3,3)
    plot(samples(:,1))
    hold on
    plot(get(gca,'Xlim'),[0.3,0.3],'r')
    ylabel('w')

    % --- alphas
    subplot(3,3,4)
    xi=0:0.005:1;
    f=ksdensity(samples(:,2),xi,'support','positive');
    plot(xi,f)
    hold on
    plot([0.03,0.03],get(gca,'Ylim'),'r')
    xlabel('alpha')

    subplot(3,3,5)
    hist(samples(:,2),500)
    hold on
    plot([0.03,0.03],get(gca,'Ylim'),'r')
    xlabel('alpha')

    subplot(3,3,6)
    plot(samples(:,2))
    hold on
    plot(get(gca,'Xlim'),[0.03,0.03],'r')
    ylabel('alpha')

    % --- t0
    subplot(3,3,7)
    xi=0:512;
    f=ksdensity(samples(:,3),xi,'support','positive');
    plot(xi,f)
    hold on
    plot([256 256],get(gca,'Ylim'),'r')
    xlabel('t0')

    subplot(3,3,8)
    hist(samples(:,3),500)
    hold on
    plot([256, 256],get(gca,'Ylim'),'r')
    xlabel('t0')

    subplot(3,3,9)
    plot(samples(:,3))
    hold on
    plot(get(gca,'Xlim'),[256 256],'r')
    ylabel('t0')
end