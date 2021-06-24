% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% This function calculates PSD in 4s windows. 
% Last updated 22.06.2021

function varargout = 6B_el_plot_spec_SSRI(X,Fs,F_max,name,varargin)
freqmark = [];

if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for j = 1:2:(length(varargin)-1)
        if ~ischar (varargin{j})
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        switch lower (varargin{j})
            case 'freqmark'
                freqmark = varargin{j+1};  
        end
    end
end
%%
[X_spec,F] = pwelch(X,1*Fs,0,4*Fs,Fs);
X_spec = X_spec(F<=F_max,:);
F = F(F<=F_max);
varargout{1} = X_spec;
varargout{2} = F;
hold on
for n = 1:size(X_spec,2)
    command = [ 'disp(''x ' num2str(n) ''')' ];
    if exist('name') == 1
        pl(n) = plot(F,10*log10(X_spec(:,n)),'linewidth',1.5, 'ButtonDownFcn', command), title(name);
    else
        pl(n) = plot(F,10*log10(X_spec(:,n)),'linewidth',1.5, 'ButtonDownFcn', command)
    end
end
yLimit_spec = get(gca,'YLim');
for k = 1:numel(freqmark)
    hold on, 
    plot(freqmark(k)*[1,1],yLimit_spec,'--','color',rand(1,3));
end

end