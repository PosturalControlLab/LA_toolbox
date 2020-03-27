function Spar =  WinterTable(m_Body, h_Body)
% ind=1; input is m_Body (kg) and h_Body (m);


        Spar.l_B = h_Body;
        Spar.m_B = m_Body;
        
        Spar.m_F = 2* 0.0145 * m_Body;
        Spar.m_L = 2* (0.0465+0.1) * m_Body; % relative mass of shank+thigh
        Spar.m_T = 0.678 * m_Body;

        Spar.l_F = 0.039 * h_Body; % vertical length...
        Spar.l_L = 0.53 * h_Body - Spar.l_F;
        Spar.l_T = (0.87-0.53) * h_Body; % heigth of glenohumeral joint above hip

        Spar.h_L = 0.553 * Spar.l_L;
        Spar.h_T = 0.626 * (0.818-0.53) * h_Body; % heigth trunk com above hip
        Spar.h_com = (Spar.h_L.*Spar.m_L + (Spar.l_L+Spar.h_T).*Spar.m_T) ./ (Spar.m_L+Spar.m_T); %com height above ankle

        % Moment of inertia for rotation about segment COM
        Spar.J_Fc = Spar.m_F .* (0.475*Spar.l_F).^2;
        Spar.J_Lc = Spar.m_L .* (0.326*Spar.l_L).^2;
        Spar.J_Tc = Spar.m_T .* (0.496*Spar.l_T).^2;

        % Moment of inertia for rotation about joints
        Spar.J_L = Spar.J_Lc + Spar.m_L.*Spar.h_L.^2;
        Spar.J_T = Spar.J_Tc + Spar.m_T.*Spar.h_T.^2;
        Spar.J_B = Spar.J_Lc + Spar.J_Tc + Spar.m_L.*Spar.h_L.^2 + Spar.m_T.*(Spar.l_L + Spar.h_T).^2;
        
        Spar.eyesAboveGround = 0.936 * h_Body;
        Spar.shoulderAboveAnkle = (0.818-0.039) * h_Body;
        Spar.shoulderToElbow = 0.186 * h_Body;
        Spar.ElbowToFinger = (0.146 + 0.108) * h_Body;
