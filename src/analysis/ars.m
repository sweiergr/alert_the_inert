function [xopt]=ars(func, LB, UB, c, rho, maxevaluations, maxlocalmin);
global currentmin 
%{
    %% ARS - accelerated random search ; 
% "On Accelerated Random Search"
% M. J. Appel, R. LaBarre, D. Radulovic (2003), in SIAM
    This file was only used for exploratory purposes and is not used in the final estimation.

%}
LB=LB(:)' ;                                  
UB=UB(:)' ;

UB_old=UB;
LB_old=LB ;

nargs=length(LB);                      %number of parameters
 

%starting values for model parameters
 
% original setup
% initial parameter vector:
% x=LB+(UB-LB).*rand(1,nargs);         

% I CHANGED this here to starting values at midpoint of bounds
x=LB+(UB-LB).*0.5;


f=feval(func,x');     %function evaluation with parameters x

xopt=x
fstar=f


r=1;
L=1;    % counts the number of local minima
i=1;

while i<=maxevaluations;

        
    UB_new= min(UB_old,xopt+r.*(UB_old-LB_old));
    LB_new= max(LB_old,xopt-r.*(UB_old-LB_old));
    

    x=LB_new+(UB_new-LB_new).*rand(1, nargs);
    f=feval(func,x');
    
    if f<fstar;

	%% NOTE: original "ars" is searching for MAXIMUM of function "func"
	%% I am obviously searching for MINIMUM !

        xopt=x;
        fstar=f;
        r=1;
        xopt;
        fstar;
        LB_new=LB;
        UB_new=UB;
    else;
        r=r/c;
        if r<rho;
            r=1;
            LB_new=LB;
            UB_new=UB;
            
            if maxlocalmin>0;
                
                localmin(L,:)=xopt;
                localminfstar(L)=fstar;

                
                % Hier lasse ich mir die 'lokalen minima' 
                % in den file "localmin.out" ausdrucken und sortieren

                diary localmin.out;     % program specific
                
                disp('-----------------------------');
                
                format long
                
                %disp('all local minima:');
                % localmin
                % disp('all local fstar:');
                % localminfstar              
                disp('all [ local minima;  fstar(rescaled) ] sorted:');
                 temp=[ localmin  (localminfstar./min(localminfstar))'];
                [tempa , tempb]=size(temp);
                tempsort=sortrows(temp,[tempb])
                currentmin=tempsort(1,1:end-1);
                format short
                
                disp('-----------------------------');
                diary off;       
                
                
                
                
                if L<maxlocalmin;
                
                    L=L+1;
                    
                    % start with new parameter vector:                
                    
                    x=LB+(UB-LB).*rand(1,nargs);                              
                    
                    f=feval(func,x');    %function evaluation with parameters x
                    xopt=x;
                    fstar=f;

                else;
                    
                    i=maxevaluations;
                    
                end;
                
            end;
                
                
        end;
    end;
    
    LB_old=LB_new;
    UB_old=UB_new;
    
    i=i+1;
    
end;
 
disp('-----------------------------');
disp('--------DONE--------');
disp('-----------------------------');                

disp('all local minima:');
[ localmin   ]

disp('all local fstar:');
localminfstar'

                

[fstar,fstar_row]=min(localminfstar);
xopt=localmin(fstar_row,:)
fstar
end