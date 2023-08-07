function figo = PlotHeatMap(x,y,v,fig,range)

    % Find the figure
    try
        if(isnumeric(fig))
            fig = figure(fig);
        end
    catch
        fig = figure();
    end
    figo = fig;

    % Interpolant
    F = scatteredInterpolant(x,y,v);
    [xq,yq] = meshgrid(-4.5:0.1:4.5);
    % F.Method = 'linear';
    F.Method = 'natural';
    vq = F(xq,yq);
    
    % Further deal with unreasonable data
    [rmax,cmax] = size(vq);
    for r = 1:rmax
        for c = 1:cmax
            if vq(r,c) < 0
                vq(r,c) = 0;
            elseif vq(r,c) > max(v)
                vq(r,c) = max(v);
            end
        end
    end
    
    % Plot
    figure(fig);
    contourf(xq,yq,vq,150,'LineColor','none');      % Color map
    colorbar;                                       % Color bar
    
    try 
        isempty(range);
    catch
        range = [min(v),max(v)];
    end
    
    if ((min(v)<range(1)) || (max(v)>range(2)) || (range(1)>range(2)))
        error('range is not set properly.');
    end
    
    ColorStepSize = length(fig.Colormap);
    ColorUpper0 = [1,1,0.5];    % Yellow
    ColorLower0 = [1,1,1];      % White
    ColorUpper0 = [1,0.5,0.5];  % Pink
    
    ColorLower = Affine(ColorLower0,ColorUpper0,(min(v)-range(1))/(range(2)-range(1)));
    ColorUpper = Affine(ColorLower0,ColorUpper0,(max(v)-range(1))/(range(2)-range(1)));
    GradRed     = linspace(ColorLower(1),ColorUpper(1),ColorStepSize)';
    GradGreen   = linspace(ColorLower(2),ColorUpper(2),ColorStepSize)';
    GradBlue    = linspace(ColorLower(3),ColorUpper(3),ColorStepSize)';
    fig.Colormap = [GradRed GradGreen GradBlue];
    
   	figure;
    mesh(xq,yq,vq);
    hold on;
    plot3(x,y,v,'.');
    
end

function m = Affine(x,y,r)
    m = x*(1-r) + y*r;
end