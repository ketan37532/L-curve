function [est_x, est_y]= error_calc_with_rssi_new(s,resolution,sat_pts,sig)
sat_pts1=sat_pts;  %%%%%received three non-collinear beacon positions (neighbors)%%%%%%%%%%%%%%%%%
N=100;
d0=1;
n=3.3;
p=3;
PL_d0_dB=55;
pt_dB=[];
RSSI1=[];

%%%%%%%%%% Formation of curves%%%%%%%%%%%%%%%%
x_1 = sat_pts1(1,1);
x_2 = sat_pts1(2,1);
x_3 = sat_pts1(3,1);
y_1 = sat_pts1(1,2);
y_2 = sat_pts1(2,2);
y_3 = sat_pts1(3,2); 

temp = [];

if(x_1==x_2)
    if(y_1==y_3)
        temp(1,:) = [x_2 y_2];
        temp(2,:) = [x_1 y_1];
        temp(3,:) = [x_3 y_3];
    elseif(y_2==y_3)
        temp(1,:) = [x_1 y_1];
        temp(2,:) = [x_2 y_2];
        temp(3,:) = [x_3 y_3];
    elseif((y_1~=y_3)&&(y_2~=y_3)&&(y_2<y_1))
        temp(1,:) = [x_2 y_2];
        temp(2,:) = [x_1 y_1];
        temp(3,:) = [x_3 y_3];
    elseif((y_1~=y_3)&&(y_2~=y_3)&&(y_1<y_2))
        temp(1,:) = [x_1 y_1];
        temp(2,:) = [x_2 y_2];
        temp(3,:) = [x_3 y_3];
    end
elseif(x_2==x_3)
    if(y_1==y_2)
        temp(1,:) = [x_3 y_3];
        temp(2,:) = [x_2 y_2];
        temp(3,:) = [x_1 y_1];
    elseif(y_3==y_1)
        temp(1,:) = [x_2 y_2];
        temp(2,:) = [x_3 y_3];
        temp(3,:) = [x_1 y_1];
    elseif((y_2~=y_1)&&(y_3~=y_1)&&(y_3<y_2))
        temp(1,:) = [x_3 y_3];
        temp(2,:) = [x_2 y_2];
        temp(3,:) = [x_1 y_1];
    elseif((y_2~=y_1)&&(y_3~=y_1)&&(y_2<y_3))
        temp(1,:) = [x_2 y_2];
        temp(2,:) = [x_3 y_3];
        temp(3,:) = [x_1 y_1];    
    end
elseif(x_1==x_3)
    if(y_1==y_2)
        temp(1,:) = [x_3 y_3];
        temp(2,:) = [x_1 y_1];
        temp(3,:) = [x_2 y_2];
    elseif(y_2==y_3)
        temp(1,:) = [x_1 y_1];
        temp(2,:) = [x_3 y_3];
        temp(3,:) = [x_2 y_2];
    elseif((y_1~=y_2)&&(y_3~=y_2)&&(y_3<y_1))
        temp(1,:) = [x_3 y_3];
        temp(2,:) = [x_1 y_1];
        temp(3,:) = [x_2 y_2];
    elseif((y_1~=y_2)&&(y_3~=y_2)&&(y_1<y_3))
        temp(1,:) = [x_1 y_1];
        temp(2,:) = [x_3 y_3];
        temp(3,:) = [x_2 y_2];    
    end
elseif((y_1==y_2)&&(x_3~=x_1)&&(x_3~=x_2)&&(x_1<x_2))
    temp(1,:) = [x_3 y_3];
    temp(2,:) = [x_1 y_1];
    temp(3,:) = [x_2 y_2];
elseif((y_1==y_2)&&(x_3~=x_1)&&(x_3~=x_2)&&(x_2<x_1))
    temp(1,:) = [x_3 y_3];
    temp(2,:) = [x_2 y_2];
    temp(3,:) = [x_1 y_1];    
elseif((y_2==y_3)&&(x_1~=x_2)&&(x_1~=x_3)&&(x_2<x_3))
    temp(1,:) = [x_1 y_1];
    temp(2,:) = [x_2 y_2];
    temp(3,:) = [x_3 y_3];
elseif((y_2==y_3)&&(x_1~=x_2)&&(x_1~=x_3)&&(x_3<x_2))
    temp(1,:) = [x_1 y_1];
    temp(2,:) = [x_3 y_3];
    temp(3,:) = [x_2 y_2];
elseif((y_1==y_3)&&(x_2~=x_1)&&(x_2~=x_3)&&(x_1<x_3))
    temp(1,:) = [x_2 y_2];
    temp(2,:) = [x_1 y_1];
    temp(3,:) = [x_3 y_3];
elseif((y_1==y_3)&&(x_2~=x_1)&&(x_2~=x_3)&&(x_3<x_1))
    temp(1,:) = [x_2 y_2];
    temp(2,:) = [x_3 y_3];
    temp(3,:) = [x_1 y_1];    
end

sat_pts1;
temp;

if(temp(1,2)>temp(2,2) && temp(2,1)<temp(3,1))
    if((temp(1,1)==temp(2,1))&&(temp(2,2)==temp(3,2)))
    curve_no = 1;
    elseif(temp(2,2)==temp(3,2))&&(temp(1,1)~=temp(2,1))&&(temp(1,1)~=temp(3,1))
    curve_no = 1.1;
    elseif(temp(1,1)==temp(2,1))&&(temp(2,2)~=temp(3,2))&&(temp(1,2)~=temp(3,2))
    curve_no = 3.2;
    end
elseif(temp(1,2)>temp(2,2) && temp(3,1)<temp(2,1))
    if(temp(1,1)==temp(2,1) && temp(3,2)==temp(2,2))
    curve_no = 2;
    elseif(temp(2,2)==temp(3,2))&&(temp(1,1)~=temp(2,1))&&(temp(1,1)~=temp(3,1))
    curve_no = 1.1;
    elseif(temp(1,1)==temp(2,1))&&(temp(2,2)~=temp(3,2))&&(temp(1,2)~=temp(3,2))
    curve_no = 4.2;
    end
elseif(temp(1,2)<temp(2,2) && temp(2,1)<temp(3,1))
    if(temp(1,1)==temp(2,1) && temp(3,2)==temp(2,2))
    curve_no = 3;
    elseif(temp(2,2)==temp(3,2))&&(temp(1,1)~=temp(2,1))&&(temp(1,1)~=temp(3,1))
    curve_no = 3.1;
    elseif(temp(1,1)==temp(2,1))&&(temp(2,2)~=temp(3,2))&&(temp(1,2)~=temp(3,2))
    curve_no = 3.2;
    end
elseif(temp(1,2)<temp(2,2) && temp(3,1)<temp(2,1))
    if(temp(1,1)==temp(2,1) && temp(3,2)==temp(2,2))
    curve_no = 4;
    elseif(temp(2,2)==temp(3,2))&&(temp(1,1)~=temp(2,1))&&(temp(1,1)~=temp(3,1))
    curve_no = 3.1;
    elseif(temp(1,1)==temp(2,1))&&(temp(2,2)~=temp(3,2))&&(temp(1,2)~=temp(3,2))
    curve_no = 4.2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%Calculation of distances based on RSSI %%%%%%%%%%%%%%%%%%%

for j=1:p
    pt_dB(j)=10-30;
end


for j=1:p
     RSSI1=[];
     d(j) = sqrt(((s(1,1)-temp(j,1))^2) +((s(1,2)-temp(j,2))^2));
     if(d(j) ==0)
        d(j) = 0.1;
     end
        for i=1:N
             sig_dB = sig;   %4
             sig_abs=sig_dB*(log(10)/10);
             x_abs=normrnd(0,sig_abs);
             x_dB=x_abs*(10/log(10));
             PL_dB(i)=PL_d0_dB+10*n*log10(d(j)/d0)+x_dB;
             A=pt_dB(j)-PL_d0_dB;
             RSSI=pt_dB(j)-PL_dB(i);
             RSSI1=[RSSI1;RSSI];
        end
        RSSI_avg=0;
        for ii=1:1:N
             RSSI_avg=RSSI_avg+RSSI1(ii);
        end
        RSSI_avg=RSSI_avg/N;
        dist(j)=10^((A-RSSI_avg)/(10*n));     
end

x_1 = temp(1,1);
x_2 = temp(2,1);
x_3 = temp(3,1);
y_1 = temp(1,2);
y_2 = temp(2,2);
y_3 = temp(3,2);

d_1 = dist(1,1);
d_2 = dist(1,2);
d_3 = dist(1,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%Algorithm 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pt_dB;
d;
dist;
sat_pts1;
est_sen=[];



d_12 = sqrt((y_1-y_2)^2+(x_1-x_2)^2);
d_23 = sqrt((y_2-y_3)^2+(x_2-x_3)^2);
d_13 = sqrt((y_1-y_3)^2+(x_1-x_3)^2);

beta_0 = acosd((d_23^2 + d_2^2 - d_3^2)/(2*d_23*d_2));
beta_1 = acosd((d_12^2 + d_2^2 - d_1^2)/(2*d_12*d_2));

tf_0 = isreal(beta_0);
tf_1 = isreal(beta_1);


%%%%%%%%%%%%%%%%%Update Distance step %%%%%%%%%%%%%%%%%%%%%%%%
if(tf_0==0)&&(tf_1==0)  
            [xout3,yout3] = circcirc(x_1,y_1,d_1,x_3,y_3,d_3);
            d_n = sqrt(((x_2-xout3(1,1))^2) +((y_2-yout3(1,1))^2));
            d_m = sqrt(((x_2-xout3(1,2))^2) +((y_2-yout3(1,2))^2));
            if ((abs(d_2-d_n)<abs(d_2-d_m)))
                num = 1;
            else 
                num = 2;
            end
            d_2 = sqrt(((x_2-xout3(1,num))^2) +((y_2-yout3(1,num))^2));

elseif(tf_0==0)  
          [xout1,yout1] = circcirc(x_1,y_1,d_1,x_2,y_2,d_2); 
          d_n = sqrt(((x_3-xout1(1,1))^2) +((y_3-yout1(1,1))^2));
          d_m = sqrt(((x_3-xout1(1,2))^2) +((y_3-yout1(1,2))^2));
          if ((abs(d_3-d_n)<abs(d_3-d_m)))
            num = 1;
          else 
            num = 2;
          end
          d_3 = sqrt(((x_3-xout1(1,num))^2) +((y_3-yout1(1,num))^2));
elseif(tf_1==0)  
          [xout2,yout2] = circcirc(x_2,y_2,d_2,x_3,y_3,d_3);
          d_n = sqrt(((x_1-xout2(1,1))^2) +((y_1-yout2(1,1))^2));
            d_m = sqrt(((x_1-xout2(1,2))^2) +((y_1-yout2(1,2))^2));
            if ((abs(d_1-d_n)<abs(d_1-d_m)))
                num = 1;
            else 
                num = 2;
            end
            d_1 = sqrt(((x_1-xout2(1,num))^2) +((y_1-yout2(1,num))^2));
end
beta_0 = acosd((d_23^2 + d_2^2 - d_3^2)/(2*d_23*d_2));
beta_1 = acosd((d_12^2 + d_2^2 - d_1^2)/(2*d_12*d_2));
alpha = acosd((d_12^2+d_23^2-d_13^2)/(2*d_12*d_23));
cos_beta_0 = cosd(beta_0);
cos_beta_1 = cosd(beta_1);
sin_beta_0 = sind(beta_0);
sin_beta_1 = sind(beta_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%Position estimation based on Curves%%%%%%%%%%%%%%%%%%%%%
    tf_0= isreal(beta_0);
    tf_1=isreal(beta_1);
    if(tf_0==0||tf_1==0)
        est_x=0;    %%%%%%%%%%if sensor can not be located%%%%%%%%%%%
        est_y=0;
    elseif(curve_no == 1)
        if(cos_beta_0 < 0 && cos_beta_1 < 0)
                x = d_2 * cosd(180-beta_0); 
                y = d_2 * sind(180-beta_0);
                est_x = x_2 - x;
                est_y = y_2 - y; 
        elseif(cos_beta_0 > 0 && cos_beta_1 < 0)
                x = d_2 * cos_beta_0;
                y = d_2 * sin_beta_0;
                est_x = x_2 + x;
                est_y = y_2 - y;
        elseif(cos_beta_0 > 0 && cos_beta_1 > 0)
                x = d_2 * cos_beta_0;
                y = d_2 * sin_beta_0;
                est_x = x_2 + x;
                est_y = y_2 + y;
        elseif(cos_beta_0 < 0 && cos_beta_1 > 0)
                x = d_2 * sin_beta_1;
                y = d_2 * cos_beta_1;
                est_x = x_2 - x;
                est_y = y_2 + y;
        elseif((cos_beta_0 == 1)) 
                est_y = y_2;
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    est_x =x_2 + d_2;
                elseif((d_3>d_23)&&(d_2<d_3))
                    est_x =x_2 - d_2;
                end
        elseif((cos_beta_1 == 1)) 
                est_x = x_2;
                if((d_2>=d_12)&&(d_1<d_2))||((d_2<d_12)&&(d_1<d_12))
                    est_y =y_2 + d_2;
                elseif((d_1>d_12)&&(d_2<d_1))
                    est_y =y_2 - d_2;
                end
        else
            est_x=0;
            est_y=0;        
        end
    elseif(curve_no== 1.1)
        if((cos_beta_0>0) && (cos_beta_1>0)) && (((beta_0<alpha) && (beta_1<alpha)) ||(beta_0>beta_1))
            x = d_2 * cos_beta_0;
            y = d_2 * sin_beta_0;
            est_x = x_2 + x;
            est_y = y_2 + y;
        elseif((cos_beta_0<0)&&(cos_beta_1>0)&&(beta_0>beta_1))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 + y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_0>beta_1)&&(beta_1<(180-alpha))&&(beta_0<180))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 + y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1>(180-alpha))&&(beta_0>(180-alpha)))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 - y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 - y;
        elseif((cos_beta_0>0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * cos_beta_0; 
            y = d_2 * sin_beta_0;
            est_x = x_2 + x;
            est_y = y_2 - y;
        elseif((cos_beta_0>0)&&(cos_beta_1>0)&&(beta_1>beta_0))
            x = d_2 * cos_beta_0; 
            y = d_2 * sin_beta_0;
            est_x = x_2 + x;
            est_y = y_2 - y;
        elseif(cos_beta_0 == 1) 
                est_y = y_2;
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    est_x =x_2 + d_2;
                elseif((d_3>d_23)&&(d_2<d_3))
                    est_x =x_2 - d_2;
                end
        elseif(cos_beta_1 == 1) 
                if((d_2>=d_12)&&(d_1<d_2))||((d_1<d_12)&&(d_2<d_12))
                    x = d_2 * cosd(beta_0);
                    y = d_2 * sind(beta_0);
                    est_x = x_2 + x;
                    est_y = y_2 + y;
                elseif((d_1>d_12)&&(d_2<d_1))
                    x = d_2 * cosd(180-beta_0);
                    y = d_2 * sind(180-beta_0);
                    est_x = x_2 - x;
                    est_y = y_2 - y;
                end
        else
            est_x=0;
            est_y=0;        
        end
    elseif(curve_no == 2)
        if(cos_beta_0 < 0 && cos_beta_1 < 0)
                x = d_2 * cosd(180-beta_0); 
                y = d_2 * sind(180-beta_0);
                est_x = x_2 + x;
                est_y = y_2 - y; 
        elseif(cos_beta_0 > 0 && cos_beta_1 < 0)
                x = d_2 * cos_beta_0;
                y = d_2 * sin_beta_0;
                est_x = x_2 - x;
                est_y = y_2 - y;
        elseif(cos_beta_0 > 0 && cos_beta_1 > 0)
                y = d_2 * sin_beta_0;
                x = d_2 * cos_beta_0;
                est_x = x_2 - x;
                est_y = y_2 + y;
        elseif(cos_beta_0 < 0 && cos_beta_1 > 0)
                x = d_2 * sin_beta_1;
                y = d_2 * cos_beta_1;
                est_x = x_2 + x;
                est_y = y_2 + y;
        elseif((cos_beta_0 == 1)) 
                est_y = y_2;
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    est_x =x_2 - d_2;
                elseif((d_3>d_23)&&(d_2<d_3))
                    est_x =x_2 + d_2;
                end
        elseif((cos_beta_1 == 1)) 
                est_x = x_2;
                if((d_2>=d_12)&&(d_1<d_2))||((d_2<d_12)&&(d_1<d_12))
                    est_y =y_2 + d_2;
                elseif((d_1>d_12)&&(d_2<d_1))
                    est_y =y_2 - d_2;
                end
        else
            est_x=0;
            est_y=0;        
        end
    elseif(curve_no == 3)
        if(cos_beta_0 < 0 && cos_beta_1 < 0)
                x = d_2 * cosd(180-beta_0); 
                y = d_2 * sind(180-beta_0);
                est_x = x_2 - x;
                est_y = y_2 + y; 
        elseif(cos_beta_0 > 0 && cos_beta_1 < 0)
                x = d_2 * cos_beta_0;
                y = d_2 * sin_beta_0;
                est_x = x_2 + x;
                est_y = y_2 + y;
        elseif(cos_beta_0 > 0 && cos_beta_1 > 0)
                y = d_2 * sin_beta_0;
                x = d_2 * cos_beta_0;
                est_x = x_2 + x;
                est_y = y_2 - y;
        elseif(cos_beta_0 < 0 && cos_beta_1 > 0)
                x = d_2 * sind(beta_1);
                y = d_2 * cosd(beta_1);
                est_x = x_2 - x;
                est_y = y_2 - y;
       elseif((cos_beta_0 == 1)) 
                est_y = y_2;
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    est_x =x_2 + d_2;
                elseif((d_3>d_23)&&(d_2<d_3))
                    est_x =x_2 - d_2;
                end
        elseif((cos_beta_1 == 1)) 
                est_x = x_2;
                if((d_2>=d_12)&&(d_1<d_2))||((d_2<d_12)&&(d_1<d_12))
                    est_y =y_2 - d_2;
                elseif((d_1>d_12)&&(d_2<d_1))
                    est_y =y_2 + d_2;
                end
        else
            est_x=0;
            est_y=0;        
        end
        
    elseif(curve_no== 3.1)
        if((cos_beta_0>0) && (cos_beta_1>0)) && (((beta_0<alpha) && (beta_1<alpha)) ||(beta_0>beta_1))
            x = d_2 * cos_beta_0;
            y = d_2 * sin_beta_0;
            est_x = x_2 + x;
            est_y = y_2 - y;
        elseif((cos_beta_0<0)&&(cos_beta_1>0)&&(beta_0>beta_1))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 - y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_0>beta_1)&&(beta_1<(180-alpha))&&(beta_0<180))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 - y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1>(180-alpha))&&(beta_0>(180-alpha)))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 + y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * cosd(180-beta_0); 
            y = d_2 * sind(180-beta_0);
            est_x = x_2 - x;
            est_y = y_2 + y;
        elseif((cos_beta_0>0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * cos_beta_0; 
            y = d_2 * sin_beta_0;
            est_x = x_2 + x;
            est_y = y_2 + y;
        elseif((cos_beta_0>0)&&(cos_beta_1>0)&&(beta_1>beta_0))
            x = d_2 * cos_beta_0; 
            y = d_2 * sin_beta_0;
            est_x = x_2 + x;
            est_y = y_2 + y;
        elseif(cos_beta_0 == 1) 
                est_y = y_2;
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    est_x =x_2 + d_2;
                elseif((d_3>d_23)&&(d_2<d_3))
                    est_x =x_2 - d_2;
                end
        elseif(cos_beta_1 == 1) 
                if((d_2>=d_12)&&(d_1<d_2))||((d_1<d_12)&&(d_2<d_12))
                    x = d_2 * cosd(beta_0);
                    y = d_2 * sind(beta_0);
                    est_x = x_2 + x;
                    est_y = y_2 - y;
                elseif((d_1>d_12)&&(d_2<d_1))
                    x = d_2 * cosd(180-beta_0);
                    y = d_2 * sind(180-beta_0);
                    est_x = x_2 - x;
                    est_y = y_2 + y;
                end
        else
            est_x=0;
            est_y=0;        
        end
        
    elseif(curve_no== 3.2)
        if((cos_beta_0>0) && (cos_beta_1>0)) && (((beta_0<alpha) && (beta_1<alpha)) ||(beta_1>beta_0))
            x = d_2 * sin_beta_1;
            y = d_2 * cos_beta_1;
            est_x = x_2 + x;
            est_y = y_2 - y;
        elseif((cos_beta_0>0)&&(cos_beta_1>0)&&(beta_0>beta_1))
            x = d_2 * sin_beta_1; 
            y = d_2 * cos_beta_1;
            est_x = x_2 - x;
            est_y = y_2 - y; 
        elseif((cos_beta_0<0)&&(cos_beta_1>0)&&(beta_0>beta_1))
            x = d_2 * sin_beta_1; 
            y = d_2 * cos_beta_1;
            est_x = x_2 - x;
            est_y = y_2 - y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1<(180-alpha))&&(beta_0>beta_1))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 - x;
            est_y = y_2 + y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_0>(180-alpha))&&(beta_1>(180-alpha))&&(beta_1<180))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 - x;
            est_y = y_2 + y;
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 + x;
            est_y = y_2 + y;
        elseif((cos_beta_0>0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 + x;
            est_y = y_2 + y;
        elseif(cos_beta_0 == 1) 
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    x = d_2 * sind(beta_1);
                    y = d_2 * cosd(beta_1);
                    est_x = x_2 + x;
                    est_y = y_2 - y;
                elseif((d_3>d_23)&&(d_2<d_3))
                    x = d_2 * sind(180-beta_1);
                    y = d_2 * cosd(180-beta_1);
                    est_x = x_2 - x;
                    est_y = y_2 + y;
                end
        elseif(cos_beta_1 == 1)
                if((d_2>=d_12)&&(d_1<d_2))||((d_1<d_12)&&(d_2<d_12))
                    est_x = x_2;
                    est_y = y_2 - d_2;
                elseif((d_1>d_12)&&(d_2<d_1))
                    est_x = x_2;
                    est_y = y_2 + d_2;
                end
        else
            est_x=0;
            est_y=0;        
        end
            
    elseif(curve_no == 4)
        if(cos_beta_0 < 0 && cos_beta_1 < 0)
                x = d_2 * cosd(180-beta_0); 
                y = d_2 * sind(180-beta_0);
                est_x = x_2 + x;
                est_y = y_2 + y; 
        elseif(cos_beta_0 > 0 && cos_beta_1 < 0)
                x = d_2 * cos_beta_0;
                y = d_2 * sin_beta_0;
                est_x = x_2 - x;
                est_y = y_2 + y;
        elseif(cos_beta_0 > 0 && cos_beta_1 > 0)
                y = d_2 * sin_beta_0;
                x = d_2 * cos_beta_0;
                est_x = x_2 - x;
                est_y = y_2 - y;
        elseif(cos_beta_0 < 0 && cos_beta_1 > 0)
                x = d_2 * sin_beta_1;
                y = d_2 * cos_beta_1;
                est_x = x_2 + x;
                est_y = y_2 - y;
        elseif((cos_beta_0 == 1)) 
                est_y = y_2;
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    est_x =x_2 - d_2;
                elseif((d_3>d_23)&&(d_2<d_3))
                    est_x =x_2 + d_2;
                end
        elseif((cos_beta_1 == 1)) 
                est_x = x_2;
                if((d_2>=d_12)&&(d_1<d_2))||((d_2<d_12)&&(d_1<d_12))
                    est_y =y_2 - d_2;
                elseif((d_1>d_12)&&(d_2<d_1))
                    est_y =y_2 + d_2;
                end
        else
            est_x=0;
            est_y=0;        
        end
     elseif(curve_no== 4.2)
         if((cos_beta_0>0) && (cos_beta_1>0)) && (((beta_0<alpha) && (beta_1<alpha)) ||(beta_1>beta_0))
            x = d_2 * sin_beta_1;
            y = d_2 * cos_beta_1;
            est_x = x_2 - x;
            est_y = y_2 - y;
        elseif((cos_beta_0>0)&&(cos_beta_1>0)&&(beta_0>beta_1))
            x = d_2 * sin_beta_1; 
            y = d_2 * cos_beta_1;
            est_x = x_2 + x;
            est_y = y_2 - y; 
        elseif((cos_beta_0<0)&&(cos_beta_1>0)&&(beta_0>beta_1))
            x = d_2 * sin_beta_1; 
            y = d_2 * cos_beta_1;
            est_x = x_2 + x;
            est_y = y_2 - y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1<(180-alpha))&&(beta_0>beta_1))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 + x;
            est_y = y_2 + y; 
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_0>(180-alpha))&&(beta_1>(180-alpha))&&(beta_1<180))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 + x;
            est_y = y_2 + y;
        elseif((cos_beta_0<0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 - x;
            est_y = y_2 + y;
        elseif((cos_beta_0>0)&&(cos_beta_1<0)&&(beta_1>beta_0))
            x = d_2 * sind(180-beta_1); 
            y = d_2 * cosd(180-beta_1);
            est_x = x_2 - x;
            est_y = y_2 + y;
        elseif(cos_beta_0 == 1) 
                if((d_2>=d_23)&&(d_3<d_2))||((d_2<d_23)&&(d_3<d_23))
                    x = d_2 * sind(beta_1);
                    y = d_2 * cosd(beta_1);
                    est_x = x_2 - x;
                    est_y = y_2 - y;
                elseif((d_3>d_23)&&(d_2<d_3))
                    x = d_2 * sind(180-beta_1);
                    y = d_2 * cosd(180-beta_1);
                    est_x = x_2 + x;
                    est_y = y_2 + y;
                end
        elseif(cos_beta_1 == 0)
                if((d_2>=d_12)&&(d_1<d_2))||((d_1<d_12)&&(d_2<d_12))
                    est_x = x_2;
                    est_y = y_2 - d_2;
                elseif((d_1>d_12)&&(d_2<d_1))
                    est_x = x_2;
                    est_y = y_2 + d_2;
                end
         else
            est_x=0;
            est_y=0;       
         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end