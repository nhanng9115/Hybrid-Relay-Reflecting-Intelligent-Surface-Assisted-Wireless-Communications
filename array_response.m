function a = array_response(azi,ele,N,array_type)

if strcmp(array_type,'UPA')
    for m= 0:sqrt(N)-1
        for n= 0:sqrt(N)-1
            a(m*(sqrt(N))+n+1) = exp( 1i* pi* ( m*sin(azi)*sin(ele) + n*cos(ele) ) );
        end
    end
else
    for n= 0:N-1
        a(n+1) = exp( 1i* pi* ( n*sin(azi) ) );
    end
end
a = a.'/sqrt(N);
end