%
% Function: updateplot(var)
%
% Updates plots on SPECGUI
%

function updateplot(pophandle,plothandle)
    val = get(pophandle,'Value');
    switch val;
        case 1        % Gfit
            Gstm = load('Gfit.dat');
            Gst  = load('Gst.dat');
            p1 =  loglog(Gst(:,1),Gst(:,2),'go',Gstm(:,1),Gstm(:,2),'k-',...
                Gst(:,1),Gst(:,3),'go',Gstm(:,1), Gstm(:,3),'k-',...
                'Linewidth',2,'Parent',plothandle);
            set(get(plothandle,'XLabel'),'String','w');
            set(get(plothandle,'YLabel'),'String','G*')
        case 2        % H
            data = load('H.dat');
            p1   = semilogx(data(:,1),data(:,2),'o-','Parent',plothandle);
            set(get(plothandle,'XLabel'),'String','s');
            set(get(plothandle,'YLabel'),'String','H(s)')
        case 3        % rho-eta
            data = load('rho-eta.dat');
            p1 =  loglog(data(:,2),data(:,3),'o-','Parent',plothandle);
            set(get(plothandle,'XLabel'),'String','rho');
            set(get(plothandle,'YLabel'),'String','eta')
            
        case 4        % Gfitd
            Gstm = load('Gfitd.dat');
            Gst  = load('Gst.dat');
            p1 =  loglog(Gst(:,1),Gst(:,2),'go',Gstm(:,1),Gstm(:,2),'k-',...
                Gst(:,1),Gst(:,3),'go',Gstm(:,1), Gstm(:,3),'k-',...
                'Linewidth',2,'Parent',plothandle);
            set(get(plothandle,'XLabel'),'String','w');
            set(get(plothandle,'YLabel'),'String','G*')

        case 5
            
            data  = load('Nopt.dat');
            Nv    = data(:,1);
            ev    = data(:,2);
            p1    = plot(Nv,ev,'o-','Parent',plothandle);
            set(get(plothandle,'XLabel'),'String','N');
            set(get(plothandle,'YLabel'),'String','error')
            
        case 6
            
            data  = load('Nopt.dat');
            Nv    = data(:,1);
            condN = data(:,3);
            p1    = plot(Nv,log10(condN),'o-','Parent',plothandle);
            set(get(plothandle,'XLabel'),'String','N');
            set(get(plothandle,'YLabel'),'String','log10(cond)')
            
        case 7
            
            data  = load('dmodes.dat');
            tau   = data(:,3);
            g     = data(:,2);
            p1    = loglog(tau,g,'o','Parent',plothandle);
            set(get(plothandle,'XLabel'),'String','tau');
            set(get(plothandle,'YLabel'),'String','g')         
    end
end