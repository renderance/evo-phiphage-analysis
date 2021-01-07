function glenn_invasibility_plotter(wortel,knol)

% Function that plots relative growths and thus invasibility for a range of phenotypes, using an internal solver function.

%{
    Performs 100 invasions onto each of the phenotypes defined by a maximal
    lysogeny propensity as given in the wortel input vector and the regulator
    threshold values given in the knol double. It does so using standard 
    parameter values. The invader phenotypes are always non-regulators with 
    propensities for lysogeny ranging from 0 to 1. Pretty nice.

    It is possible to perform an invasion onto a non-regulator by setting
    knol to zero. 

    My apologies for the nonsensical comments, but it was one of those
    days...
%}

    % Prepare some parameters:
    bi=1;
    de=.3;
    al=.4;
    be=.25;
    ga=4;
    dp=.5;
    t = 0:.5:200;


    for jj = wortel
        % Prepare 100 invader phenotypes, stick them in a vector called 'a'
        a=linspace(0,1,100);
        
        % Open a figure, you know, just in case you'll need it.
        figure;

        % Latex is nice, so if you ever need to put text in the figure,
        % make it be interpreted as latex text. Just because you can.
        set(0,'defaultTextInterpreter','latex');

        % Remember 'a'? Yeah, we're gonna iterate over its entries.
        for ii = 1:length(a)
            
            % Let's make an RGB color using a range like the iterator.
            color = [1-(ii-1)/length(a) 0 (ii-1)/length(a)];

            % Calm the user with some text, in case he gets worried you've
            % been wasting his time with useless color mixing and window opening. 
            % You have, but (s)he doesn't need to know that.
            fprintf(['Working on invasion number ',num2str(ii,'%i'),':'])

            % Create a vector with the current resident jj and invader
            % a(ii). I blame my mother for my creative variable names.
            ps=[jj a(ii)];

            % Create starting condition vector.
            init_x=[0.7; 0; 0; 1e-2; 1e-2];

            % Solve the five ordinary differential equations, using the
            % invasibility solver, because that neglects the influence of
            % the invader on the sensitive population. Way epic.
            y = glenn_solver_invasibility(bi,de,al,be,ga,dp,ps,t,init_x);

            % Calculate the number of resident virions and lysogens.
            z = y(:,2)+y(:,4);
            % Calculate the number of invader virions and lysogens.
            w = y(:,3)+y(:,5);

            % Steal a penguin from the zoo and feed it a graphical object.
            penguin=subplot(2,1,2);
            % Put invader density relative to resident in that graphical
            % object that you just stuffed your penguin with.
            plot(t,w./z,'Color',color,'Parent',penguin)
            % Now hold on to your penguin, for it will try to flee.
            hold(penguin,'on')

            % Your penguin is a humbolt, so it needs another graphical
            % object. Stuff it in there.
            humbolt=subplot(2,1,1);
            % That graphical object needs some more entries. Duplicate the
            % figure from before. You'll need that so you can zoom in on a
            % part of it, later.
            plot(t,w./z,'Color',color,'Parent',humbolt)
            % Oh, yeah. Don't let go of it, or it will forget the data.
            % Keep that swimming bird from forgetting your precious plots.
            hold(humbolt,'on')

        end

        % Add some horizontal lines to your penguins.
        plot(t,z./z,'Color','g','Parent',humbolt);
        plot(t,z./z,'Color','g','Parent',penguin);
        % Zoom a part of your penguin into a particular window of data.
        xlim(penguin,[90 130]);

        % Make sure everybody knows your penguins are yours, so label them.
        ylabel(humbolt,"$\frac{P_{inv}+J_{inv}}{P_{res}+J_{res}}$");
        xlabel(humbolt,"$Time (b^{-1})$");
        ylabel(penguin,"$\frac{P_{inv}+J_{inv}}{P_{res}+J_{res}}$");
        xlabel(penguin,"$Time (b^{-1})$");


    end




    function x = glenn_solver_invasibility(b,d,alpha,beta,gamma,delta,phis,time,x0)

    % Invoke the ode45 solver for the system of equations, and if it fails, try ode23s:
    lastwarn('');
    [t,x] = ode45(@RHS,time,x0);
    [~,msgid]=lastwarn;
        if msgid ~= 0
            [t,x] = ode23s(@RHS,time,x0);
            fprintf(' finished with ODE23S.\n')
        else
            fprintf(' finished with ODE45.\n')
        end



        % Construct the Right-Hand-Side vector of equations used by the solver:
        function fx = RHS(t,x)

            x(x<0)=0;

            % Create an empty column-vector for all the equations:
            fx = zeros(5,1);

            % Define the regulation threshold value:
            ts = knol;
            
            % Create the first index, which always contains the sensitives:
            fx(1) =                         b*x(1)*(1-x(1)-x(2))...
                                            -d*x(1)...
                                            -beta*x(1)*x(4);

            fx(2) =                         b*x(2)*(1-x(1)-x(2))...
                                            -d*x(2)...
                                            -alpha*x(2)...
                                            +heaviside(beta*x(1)*sum(x(4))-ts)*phis(1)*beta*x(1)*x(4);

            fx(4) =                         gamma*(...
                                            alpha*x(2)+(1-heaviside(beta*x(1)*sum(x(4))-ts)*phis(1))*beta*x(1)*x(4))...
                                            -delta*x(4);

            fx(3) =                         b*x(3)*(1-x(1)-x(2))...
                                            -d*x(3)...
                                            -alpha*x(3)...
                                            +phis(2)*beta*x(1)*x(5);

            fx(5) =                         gamma*(...
                                            alpha*x(3)+(1-phis(2))*beta*x(1)*x(5))...
                                            -delta*x(5);

        end



    end


end
