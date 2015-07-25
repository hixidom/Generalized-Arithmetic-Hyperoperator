function [z]=hyperop(x,n,a,cutoff) %xth unary hyperoperator applied n times on a (expansion truncated at cutoff)

    %For relatively trivial cases of both x and n being integers...
    if rem(x,1)==0 && rem(n,1)==0
        z=a;
        if n>=0
            for i=1:n
                switch x
                    case 1
                        z=z*2;
                    case 2
                        z=z^2;
                    case 3
                        z=z^z;
                    otherwise z=hyperop(x-1,log(z)/log(2),z,cutoff);

                end
            end
        %Inverse operations...
        else
            for i=1:-n
                switch x
                    case 1
                        z=z/2;
                    case 2
                        z=sqrt(z);
                    case 3
                        z=log(z)/lambertw(log(z));
                    otherwise z=hyperop(x-1,-log(z)/log(2),z,cutoff); %Probably it is not correct to simply throw in a minus sign in front of n
                end
            end
        end
        
    %For slightly less trivial case of x alone being an integer...
    elseif rem(x,1)==0
        z=0;
        for i=0:cutoff
            prod=1;
            for j=1:i
                prod=prod*(n-1+j);
            end
            sum2=0;
            for j=0:i
                sum2=sum2+((-1)^j*nchoosek(i,j)*hyperop(x,-j,a,cutoff));
            end
            z=z+(prod*sum2/factorial(i));
        end
    
    %For the completely general case of x, n, and a being real numbers...
    else
        z=0;
        for i=1:cutoff
            prod=1;
            for j=1:i
                prod=prod*(x+1-j);
            end
            sum2=0;
            for j=1:i
                sum2=sum2+((-1)^j*nchoosek(i,j)*hyperop(j,1,a,cutoff));
            end
            z=z+((-1)^i*prod*sum2/factorial(i));
        end
    end
end