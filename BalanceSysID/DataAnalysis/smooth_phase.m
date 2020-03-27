function xin = smooth_phase(xin)

    for k=2:length(xin)
        while xin(k)>xin(k-1)+180
            xin(k)=xin(k)-360;
        end
        while xin(k)<xin(k-1)-180
            xin(k)=xin(k)+360;
        end
    end

