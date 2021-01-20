
    x1 = [1 2 3 4 5 5.5 7 8 9 9.5 10];
    v1 = [0 0 0 0.5 0.4 1.2 1.2 0.1 0 0.3 0.6];
    xq1 = linspace(1,10, 1000);
    vqm = makima(x1,v1,xq1); % same as vq = interp1(x,v,xq,'makima')
    vqm2  = load('../makima_F90/data')

    figure 
    plot(x1, v1, '*k', xq1, vqm, vqm2(:,1), vqm2(:,2), ...
         'LineWidth', 3) 
