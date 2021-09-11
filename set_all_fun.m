classdef set_all_fun
    %SET_ALL_FUN 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        f1_1;
        f1_2;
        f2_1;
        f2_2;
        f3_1;
        f3_2;
        
        ad_a_x;
        ad_a_y;
        ad_m;
        as;
        B;
        argument;
        
        Ax;Ay;
        bx;by;
    end
    methods
        function obj = set_all_fun(apha,c_x,c_y,glass,m_in,m1_in,s_in,s1_in)
            %SET_FUN 构造此类的实例
            %   此处显示详细说明
            fd = 33;
            fs = 1200;
            obj_dis = 1200;
            pix_gauge = 5.5e-3;
            %垂轴放大率
            l = pix_gauge*sqrt((glass(1,1)-glass(1,2))^2 + (glass(2,1)-glass(2,2))^2);
            L = 500;%实际尺寸
            B = l/L;%垂轴放大率0.0282662
            %B = 1;%test
            obj.B = B;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ad = 1.6222e-4;%物方视场角变化1s后所成像的移动量（形变）
            as = 5.8833e-3;%物方视场角变化1s后所成像的移动量（太阳）
            %标定镜头截距
            Ld = 0.584/tan(pi/180);%33.4573    0.584是移动1°时形变镜光斑位移量
            Ls = 0.353/tan(pi/60/180);%1213.5246    0.353是移动1′时形变镜光斑位移量
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %专利方法标定中心光斑截距
            theta_x_basic = 3.2 *pi/3600/180;
            theta_y_basic = -2.2 *pi/3600/180;
            
            theta_x_1 = 86.8 *pi/3600/180;
            theta_y_1 = -0.8 *pi/3600/180;
            theta_x_2 = -82.8 *pi/3600/180;
            theta_y_2 = -0.4 *pi/3600/180;
            theta_x_3 = -212.4 *pi/3600/180;
            theta_y_3 = -0.6 *pi/3600/180;
            
            m_sam_1 = [1039.36, 1024.30];
            m_sam_2 = [1036.63, 1027.51];
            m_sam_3 = [1033.7, 1029.84];
            
            s_sam_1 = [1057.0, 936.8];
            s_sam_2 = [927.7, 1064.4];
            s_sam_3 = [829.5, 1162.5];
            
            %{
            Ld_x_1 = (m_sam_1(1) - m_in(1)) * pix_gauge / tan(theta_x_1 - theta_x_basic);
            Ld_x_2 = (m_sam_2(1) - m_in(1)) * pix_gauge / tan(theta_x_2 - theta_x_basic);
            Ld_x_3 = (m_sam_3(1) - m_in(1)) * pix_gauge / tan(theta_x_3 - theta_x_basic);
            
            Ld_y_1 = (m_sam_1(2) - m_in(2)) * pix_gauge / tan(theta_y_1 - theta_y_basic);
            Ld_y_2 = (m_sam_2(2) - m_in(2)) * pix_gauge / tan(theta_y_2 - theta_y_basic);
            Ld_y_3 = (m_sam_3(2) - m_in(2)) * pix_gauge / tan(theta_y_3 - theta_y_basic);
            
            Ls_x_1 = (s_sam_1(1) - s_in(1)) * pix_gauge / tan(theta_x_1 - theta_x_basic);
            Ls_x_2 = (s_sam_2(1) - s_in(1)) * pix_gauge / tan(theta_x_2 - theta_x_basic);
            Ls_x_3 = (s_sam_3(1) - s_in(1)) * pix_gauge / tan(theta_x_3 - theta_x_basic);
            
            Ls_y_1 = (s_sam_1(2) - s_in(2)) * pix_gauge / tan(theta_y_1 - theta_y_basic);
            Ls_y_2 = (s_sam_2(2) - s_in(2)) * pix_gauge / tan(theta_y_2 - theta_y_basic);
            Ls_y_3 = (s_sam_3(2) - s_in(2)) * pix_gauge / tan(theta_y_3 - theta_y_basic);
            
            %Ld = obj.average(Ld_x_1, Ld_x_2, Ld_x_3, Ld_y_1, Ld_y_2, Ld_y_3);
            %Ls = obj.average(Ls_x_1, Ls_x_2, Ls_x_3, Ls_y_1, Ls_y_2, Ls_y_3);
            
            %}
            
            mode = 2;
            switch mode
                case 1
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %采用垂直入射光斑作为输入值
                    Ld_1 = (sqrt( (m_sam_1(1) - m_in(1))^2+(m_sam_1(2) - m_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_1 - theta_x_basic))^2+(tan(theta_y_1 - theta_y_basic))^2 ));
                    Ld_2 = (sqrt( (m_sam_2(1) - m_in(1))^2+(m_sam_2(2) - m_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_2 - theta_x_basic))^2+(tan(theta_y_2 - theta_y_basic))^2 ));
                    Ld_3 = (sqrt( (m_sam_3(1) - m_in(1))^2+(m_sam_3(2) - m_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_3 - theta_x_basic))^2+(tan(theta_y_3 - theta_y_basic))^2 ));
                    
                    Ls_1 = (sqrt( (s_sam_1(1) - s_in(1))^2+(s_sam_1(2) - s_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_1 - theta_x_basic))^2+(tan(theta_y_1 - theta_y_basic))^2 ));
                    Ls_2 = (sqrt( (s_sam_2(1) - s_in(1))^2+(s_sam_2(2) - s_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_2 - theta_x_basic))^2+(tan(theta_y_2 - theta_y_basic))^2 ));
                    Ls_3 = (sqrt( (s_sam_3(1) - s_in(1))^2+(s_sam_3(2) - s_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_3 - theta_x_basic))^2+(tan(theta_y_3 - theta_y_basic))^2 ));
                    
                    Ld = obj.average(Ld_1, Ld_2, Ld_3);
                    Ls = obj.average(Ls_1, Ls_2, Ls_3);
                
                case 2
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %不采用标定垂直入射光斑作为输入值
                    
                    %Ld_1 = (sqrt( (m_sam_1(1) - m_in(1))^2+(m_sam_1(2) - m_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_1 - theta_x_basic))^2+(tan(theta_y_1 - theta_y_basic))^2 ));
                    Ld_2 = (sqrt( (m_sam_2(1) - m_sam_1(1))^2+(m_sam_2(2) - m_sam_1(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_2 - theta_x_1))^2+(tan(theta_y_2 - theta_y_1))^2 ));
                    Ld_3 = (sqrt( (m_sam_3(1) - m_sam_1(1))^2+(m_sam_3(2) - m_sam_1(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_3 - theta_x_1))^2+(tan(theta_y_3 - theta_y_1))^2 ));
                    
                    %Ls_1 = (sqrt( (s_sam_1(1) - s_in(1))^2+(s_sam_1(2) - s_in(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_1 - theta_x_basic))^2+(tan(theta_y_1 - theta_y_basic))^2 ));
                    Ls_2 = (sqrt( (s_sam_2(1) - s_sam_1(1))^2+(s_sam_2(2) - s_sam_1(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_2 - theta_x_1))^2+(tan(theta_y_2 - theta_y_1))^2 ));
                    Ls_3 = (sqrt( (s_sam_3(1) - s_sam_1(1))^2+(s_sam_3(2) - s_sam_1(2))^2 )) * pix_gauge / (sqrt( (tan(theta_x_3 - theta_x_1))^2+(tan(theta_y_3 - theta_y_1))^2 ));
                    
                    Ld = obj.average(Ld_2, Ld_3);
                    Ls = obj.average(Ls_2, Ls_3);

            end
            
            %%%专利定义
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %中心光斑
            Ld = 33.38;%检测中心测试主截距
            %Ld = 1200;%test
            
            ad_m = Ld;
            obj.ad_m = ad_m;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %毛玻璃
            ad_a_x_1 = Ld*(sec(atan((glass(1,1)-m_in(1))*pix_gauge/Ld)))^2;
            ad_a_x_2 = Ld*(sec(atan((glass(1,2)-m_in(1))*pix_gauge/Ld)))^2;
            ad_a_x_3 = Ld*(sec(atan((glass(1,3)-m_in(1))*pix_gauge/Ld)))^2;
            ad_a_x = obj.average(ad_a_x_1,ad_a_x_2,ad_a_x_3);    %34.1244

            obj.ad_a_x = ad_a_x;
            
            ad_a_y_1 = Ld*(sec(atan((glass(2,1)-m_in(2))*pix_gauge/Ld)))^2;
            ad_a_y_2 = Ld*(sec(atan((glass(2,2)-m_in(2))*pix_gauge/Ld)))^2;
            ad_a_y_3 = Ld*(sec(atan((glass(2,3)-m_in(2))*pix_gauge/Ld)))^2;
            ad_a_y = obj.average(ad_a_y_1,ad_a_y_2,ad_a_y_3);    %34.1115
            
            obj.ad_a_y = ad_a_y;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %太阳光斑
            as = Ls;
            obj.as = as;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m = m_in;
            m1 = m1_in;
            s = s_in;
            s1 = s1_in;
            %eps_in = 1e-4;
            x0 = obj.aver(glass(1,1),glass(1,2));
            y0 = obj.aver(glass(2,1),glass(2,2));
            sunglass_pos = obj.sunglass_pos(x0,y0,glass);
            sunglass_x = sunglass_pos(1);
            sunglass_y = sunglass_pos(2);
            syms dx dy beta_x beta_y gamma_x gamma_y;
            
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %毛玻璃光斑的移动量方程
            obj.f1_1 = (dx*B + beta_x*ad_a_x) /pix_gauge - (c_x);%像元数为单位 dx,dy为像元数
            obj.f1_2 = (dy*B + beta_y*ad_a_y) /pix_gauge - (c_y);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %形变太阳光斑的移动量方程
            obj.f2_1 = (beta_x + gamma_x)*ad_m/pix_gauge - (m1(1)-m(1));%像距为单位
            obj.f2_2 = (beta_y + gamma_y)*ad_m/pix_gauge - (m1(2)-m(2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %太阳指向光斑的移动量方程
            Tx = (sunglass_x*(cos(apha)-1) - sunglass_y*sin(apha) + dx*B/pix_gauge) * (pix_gauge / B);
            Ty = (sunglass_x*sin(apha) + sunglass_y*(cos(apha)-1) + dy*B/pix_gauge) * (pix_gauge / B);
            t=(sunglass_x*(cos(0)-1) - sunglass_y*sin(0) + dx*B/pix_gauge) * (pix_gauge / B);
            %obj.f3_1 = (gamma_x*as - Tx)/pix_gauge - (s1(1)-s(1));%注意光学倒像与实际镜头移动方向之间的关系相反，因此-Tx
            obj.f3_1 = (gamma_x*as - Tx) - (s1(1)-s(1))*pix_gauge;
            obj.f3_2 = (gamma_y*as - Ty) - (s1(2)-s(2))*pix_gauge;
            %obj.f3_1 = (sunglass_x*(cos(apha)-1) - sunglass_y*sin(apha) + dx*fd/fs/pix_gauge)*obj_dis/fd + gamma_x*fs/pix_gauge - (s1(1)-s(1));
            %obj.f3_2 = (sunglass_x*sin(apha) + sunglass_y*(cos(apha)-1) + dy*fd/fs/pix_gauge)*obj_dis/fd + gamma_y*fs/pix_gauge - (s1(2)-s(2));
            %}
            
            
            %%%%%%%%%%%%%%modify210812%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %毛玻璃光斑的移动量方程
            obj.f1_1 = (-dx*B + beta_x*ad_a_x) /pix_gauge - (c_x);%像元数为单位 dx,dy为像元数，根据倒像关系，dx取负号
            obj.f1_2 = (-dy*B + beta_y*ad_a_y) /pix_gauge - (c_y);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %形变太阳光斑的移动量方程
            obj.f2_1 = (beta_x + gamma_x)*ad_m/pix_gauge - (m1(1)-m(1));%像距为单位
            obj.f2_2 = (beta_y + gamma_y)*ad_m/pix_gauge - (m1(2)-m(2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %太阳指向光斑的移动量方程
            Tx = ( -(sunglass_x*(cos(apha)-1) - sunglass_y*sin(apha)) + dx*B/pix_gauge) * (pix_gauge / B);%-apha+dx
            Ty = ( -(sunglass_x*sin(apha) + sunglass_y*(cos(apha)-1)) + dy*B/pix_gauge) * (pix_gauge / B);
            obj.f3_1 = (gamma_x*as + Tx) - (s1(1)-s(1))*pix_gauge;
            obj.f3_2 = (gamma_y*as + Ty) - (s1(2)-s(2))*pix_gauge;
            
            
            %%%%%%out_argument%%%%%%%%%%
            obj.argument = [B, ad_a_x, ad_a_y, ad_m, as, sunglass_x, sunglass_y];
            
            %%%%%%%set_argu_matrix%%%%%%%%%%%%
            obj.Ax = [-B, ad_a_x, 0;
                0, ad_m, ad_m;
                -1, 0, as];
            obj.Ay = [-B, ad_a_y, 0;
                0, ad_m, ad_m;
                -1, 0, as];
            obj.bx = [c_x*pix_gauge;
                (m1(1)-m(1))*pix_gauge;
                (s1(1)-s(1))*pix_gauge - (sunglass_x*(cos(apha)-1) - sunglass_y*sin(apha))*(pix_gauge / B) ];
            obj.by = [c_y*pix_gauge;
                (m1(2)-m(2))*pix_gauge;
                (s1(2)-s(2))*pix_gauge - (sunglass_x*sin(apha) + sunglass_y*(cos(apha)-1))*(pix_gauge / B) ];
        end
        function f = f_out(obj)
            f = [obj.f1_1 obj.f1_2 obj.f2_1 obj.f2_2 obj.f3_1 obj.f3_2];
        end
        function argument = argument_set(obj)
            argument = [obj.B, obj.ad_x, obj.ad_y, obj.ad_m, obj.as, ojb.Tx, obj.Ty];
        end
        
        function sunglass_pos = sunglass_pos(~,x0,y0,glass)
            x = x0 + (glass(1,1) - glass(1,2))/10;%11个光栅
            y = y0 + (glass(2,1) - glass(2,2))/10;%11个光栅
            sunglass_pos = [x,y];
        end
        %{
        function sunglass_pos = cal_sunglass(obj,glass,init_sunglass_input,eps_in)
            sunglass_pos = obj.newton_sunglass(glass,init_sunglass_input,eps_in);
        end
        function cal_sunglass_output = newton_sunglass(obj,glass,init_sunglass_input,eps_in)
             %output:'sunglass_x' 'sunglass_y'
            i_value = init_sunglass_input;
            x0 = obj.aver(glass(1,1),glass(1,2));
            y0 = obj.aver(glass(2,1),glass(2,2));
            k = (glass(2,1)-glass(2,2))/(glass(1,1)-glass(1,2));
            sunglass_fun = set_sunglass_fun(x0,y0,k);
            sunglass_dfun = set_sunglass_dfun();
            for i = 1 : 1e2
                f = double( subs(sunglass_fun.f_out(),{'sunglass_x' 'sunglass_y'},{i_value(1) i_value(2)}) );
                df = double( subs(sunglass_dfun.df_out(x0,y0,k),{'sunglass_x' 'sunglass_y'},{i_value(1) i_value(2)}) );
                cal_value = i_value - f/df;
                if( abs(cal_value - i_value) < eps_in )
                    break;
                end
                i_value = cal_value; 
            end
            cal_sunglass_output = i_value;
        end
        %}
        function av = aver(~,a,b)
            av = (a+b)/2;
        end
        function av = average(~,varargin)
            a = 0;
            for i = 1:nargin-1  %包含有一个默认变量obj，未显示
                a = a + varargin{i};
            end
            av = a/(nargin-1);
        end
        
    end
end

