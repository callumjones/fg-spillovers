function VARirplot(IR,VARopt,INF,SUP)
% =======================================================================
% Plot the IRs computed with VARir
% =======================================================================
% VARirplot(IR,VARopt,vnames,INF,SUP)
% -----------------------------------------------------------------------
% INPUT
%   - IR(:,:,:) : matrix with IRF (H horizons, N variables, N shocks)
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%   - INF: lower error band
%   - SUP: upper error band
% -----------------------------------------------------------------------
% EXAMPLE
%   - See VARToolbox_Code.m in "../Primer/"
% =======================================================================
% VAR Toolbox 3.0
% Ambrogio Cesa-Bianchi
% ambrogiocesabianchi@gmail.com
% March 2012. Updated November 2020
% -----------------------------------------------------------------------


%% Check inputs
%================================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
% If there is VARopt get the vnames
vnames = VARopt.vnames;
% Check they are not empty
if isempty(vnames)
    error('You need to add label for endogenous variables in VARopt');
end

%% Retrieve and initialize variables 
%================================================
filename = [VARopt.figname 'IR_'];
quality = VARopt.quality;
suptitle = VARopt.suptitle;
pick = VARopt.pick;

% Initialize IR matrix
[nsteps, nvars, nshocks] = size(IR);

% If one shock is chosen, set the right value for nshocks
if pick<0 || pick>nvars
    error('The selected shock is non valid')
else
    if pick==0
        pick=1;
    else
        nshocks = pick;
    end
end

% Define the rows and columns for the subplots
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);

vars_to_plot = VARopt.vars_to_plot ;
sh_to_plot = VARopt.sh_to_plot ;

%% Plot
%================================================
SwatheOpt = PlotSwatheOption;
SwatheOpt.marker = '*';
SwatheOpt.trans = 1;
%FigSize(VARopt.FigSize(1),VARopt.FigSize(2))

counter = 1 ;

for jj=sh_to_plot
    for ii=vars_to_plot
        subplot(nshocks,nvars,counter);
        set(gcf, 'DefaultAxesLineWidth', 1);
        set(gcf, 'DefaultLineLineWidth', 1);
        set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
        set(gcf, 'DefaultLegendInterpreter','latex');
        set(gcf, 'DefaultAxesFontSize',6);
        %set(gcf,'PaperPosition',[0 0 11 6]) ;
        plot(steps,IR(:,ii,jj),'LineStyle','-','Color','k'); hold on
        ax = gca ;
        ax.YGrid='on' ;
        ax.MinorGridLineStyle='-' ;
        box off ;
        plot(x_axis,'-k','LineWidth',0.5); hold on
        if exist('INF','var') && exist('SUP','var')
            PlotSwathe(IR(:,ii,jj),[INF(:,ii,jj) SUP(:,ii,jj)],SwatheOpt); hold on;
        end
        xlim([1 nsteps]);
        title([vnames{ii} ' to ' vnames{jj}],'Interpreter','latex','FontSize',7.5); 
        %set(gca, 'Layer', 'top');

        counter = counter + 1 ;
    end
end
    % Save
    FigName = VARopt.FileName;
    if quality 
        if suptitle==1
            Alphabet = char('a'+(1:nshocks)-1);
            SupTitle([Alphabet(jj) ') IR to a shock to '  vnames{jj}])
        end
        set(gcf, 'Color', 'w');
        export_fig(FigName,'-pdf','-painters')
    else

        print('-depsc',FigName);
    end

close all
