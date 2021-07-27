function balancedAmplification(I,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'I')
addParameter(Parser,'tau',100)
addParameter(Parser,'W',4*(2/7)*[1 -1.1; 1 -1.1])
addParameter(Parser,'r0',[0;0])
addParameter(Parser,'dr0',[0;0])

parse(Parser,I,varargin{:})

I = Parser.Results.I;
tau = Parser.Results.tau;
W = Parser.Results.W;
r0 = Parser.Results.r0;
dr0 = Parser.Results.dr0;

%% Find response to input over time
dr = nan(size(I));
r = nan(size(I));
for ti = 1:size(I,2)
    if ti == 1
        dr(:,1) = dr0;
        r(:,1) = r0 + dr(:,1)/tau;
    else
        dr(:,ti) = -r(:,ti-1) + W*r(:,ti-1) + I(:,ti);
        r(:,ti) = r(:,ti-1) + dr(:,ti)/tau;
    end
end

%% Plot results
figure
t = (0:size(I,2)-1)/tau;
subplot(3,1,1)
plot(t,I')
subplot(3,1,2)
plot(t,r')
subplot(3,1,3)
plot(t,dr')
