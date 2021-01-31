function ans = inOpenInterval(x,a,b)
    if a>b 
        if x>b && x<a
            ans = true; 
        else 
            ans = false;
        end
    else
        if x>a && x<b
            ans = true;
        else
            ans = false;
        end
    end
end