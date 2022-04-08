class SFDandBMD():
    '''
    READ:
    1. The goal of the function is to generate Shear Force Diagram (SFD) and Bending Moment Diagram (BMD).
       The given inputs include the length of the simply supported beam, forces and locations, moment and 
       locations (if any). Force Balance and Moment Balance are used to solve the beam with given inputs.
       SFD and BMD are then plotted. IMPORTANT:The function is applicable for only for determinant simply 
       supported beams. The function cannot be applied to any statically undeterminant beam.
    2. The function is only tested by the author with limited cases and the correctness 
       of the function cannot be guaranteed. Engineering discretion needs to be used when 
       apply to practical work. The author shall not be held liable for any consequences 
       in the use of this function.
    3. Since the function is an individual work with limited coding abilities, any corrections
       or modifications suggestions are appreciated. Thank you in advance.
    '''

    def inputPara(self):
        L_mag = float(input('Enter the length of the simply supported beam (m): '))
        if L_mag <= 0:
            print('Length of the beam cannot be less or equal to 0.')
            self.inputPara()
        return L_mag
    
    def inputSupport(self, L):
        d_S1 = float(input('Enter the first support location, with left fixed end being 0 (m): '))
        if d_S1 > L or d_S1 < 0:
            print('Support is located outside of the beam range.')
            Flag = False
            while Flag == False:
                d_S1 = float(input('Enter the first support location, with left fixed end being 0 (m): '))
                if d_S1 > L or d_S1 < 0:
                    print('Support is located outside of the beam range.')
                else:
                    Flag = True
        d_S2 = float(input('Enter the second support location, with left fixed end being 0 (m): '))
        if d_S2 > L or d_S2 < 0 or d_S2 == d_S1:
            print('Support is located outside of the beam range, or has coincided with the other support.')
            Flag = False
            while Flag == False:
                d_S2 = float(input('Enter the second support location, with left fixed end being 0 (m): '))
                if d_S2 > L or d_S2 < 0 or d_S2 == d_S1:
                    print('Support is located outside of the beam range, or has coincided with the other support.')
                else:
                    Flag = True
        
        return [d_S1,d_S2]
        
    def inputForce(self, L):
        import numpy as np
        import numpy.linalg
        import sympy as sp

        X = sp.symbols('X')
        P = sp.symbols('P')
        Q = sp.symbols('Q')

        F_type = input('Enter the force type, P (Point Load) or D (Distributed Load): ')

        if F_type == 'P':
            F_mag = float(input('Enter the magnitude of the load, with negative indicating downward (N): '))
            d_F = float(input('Enter the location of the load, with left fixed end being 0 (m): '))

            if d_F > L or d_F < 0:
                print('Point force is located outside of the beam range.')
                Option = bool(int(input('Is the force still wanted, 1 (Yes) or 0 (No): ')))
                return 'False', 'False', Option, 'False', 'False', 'False', 'False'

            F_mag_oi = [F_mag]
            d_F_oi = [d_F]

            fxf = F_mag
            mxf = F_mag*X + Q

        elif F_type == 'D':
            # F_Dstbt_type = input('Enter the distributed load type, C (Constant Load), L (Linearly Variant Load), H (Linear Increase then Linear Decrease), or Q (Qudraticly Variant Load): ')
            F_Dstbt_type = input('Enter the distributed load type, C (Constant Load), L (Linearly Variant Load), or Q (Qudraticly Variant Load): ')
            if F_Dstbt_type == 'C':
                F_mag = float(input('Enter the magnitude of the constant distributed load, with negative indicating downward (N/m): '))
                d_F1 = float(input('Enter the starting location of the load, with left fixed end being 0 (m): '))
                d_F2 = float(input('Enter the ending location of the load, with left fixed end being 0 (m): '))
                
                if d_F1 > L or d_F1 < 0 or d_F2 > L or d_F2 < 0:
                    print('One or more end(s) of the distributed load is located outside of the beam range.')
                    Option = bool(int(input('Is the force still wanted, 1 (Yes) or 0 (No): ')))
                    return 'False', 'False', Option, 'False', 'False', 'False', 'False'

                fxf = F_mag*X + P
                mxf = 0.5*F_mag*X**2 + P*X + Q

                F_mag_oi = [F_mag, F_mag]
                d_F_oi = [d_F1, d_F2]

                d_F = (d_F1 + d_F2)/2
                F_mag = F_mag*(d_F2 - d_F1)

            elif F_Dstbt_type == 'L':
                F_mag1 = float(input('Enter the magnitude of the distributed load on the left end, with negative indicating downward (N/m): '))
                F_mag2 = float(input('Enter the magnitude of the distributed load on the right end, with negative indicating downward (N/m): '))
                d_F1 = float(input('Enter the starting location of the load, with left fixed end being 0 (m): '))
                d_F2 = float(input('Enter the ending location of the load, with left fixed end being 0 (m): '))

                if d_F1 > L or d_F1 < 0 or d_F2 > L or d_F2 < 0:
                    print('One or more end(s) of the distributed load is located outside of the beam range.')
                    Option = bool(int(input('Is the force still wanted, 1 (Yes) or 0 (No): ')))
                    return 'False', 'False', Option, 'False', 'False', 'False', 'False'
                
                F_mag_oi = [F_mag1, F_mag2]
                d_F_oi = [d_F1, d_F2]
                
                A = np.array([[d_F1, 1],[d_F2, 1]])
                B = np.array([[F_mag1],[F_mag2]])
                coef = np.linalg.inv(A)@B
                a, b = coef[0][0], coef[1][0]

                F_mag = (F_mag1 + F_mag2)*(d_F2 - d_F1)*0.5
                d_F = ((1/3)*a*(d_F2**3 - d_F1**3) + (1/2)*b*(d_F2**2 - d_F1**2))/F_mag

                fxf = (1/2)*a*X**2 + b*X + P
                mxf = (1/3)*(1/2)*a*X**3 + (1/2)*b*X**2 + P*X + Q

            elif F_Dstbt_type == 'Q':
                F_mag1 = float(input('Enter the magnitude of the distributed load on the first location, with negative indicating downward (N/m): '))
                F_mag2 = float(input('Enter the magnitude of the distributed load of the second location, with negative indicating downward (N/m): '))
                F_mag3 = float(input('Enter the magnitude of the distributed load on the last location, with negative indicating downward (N/m): '))
                
                d_F1 = float(input('Enter the first location of the load, with left fixed end being 0 (m): '))
                d_F2 = float(input('Enter the second location of the load, with left fixed end being 0 (m): '))
                d_F3 = float(input('Enter the last location of the load, with left fixed end being 0 (m): '))

                if d_F1 > L or d_F1 < 0 or d_F2 > L or d_F2 < 0 or d_F3 > L or d_F3 < 0:
                    print('One or more end(s) of the distributed load is located outside of the beam range.')
                    Option = bool(int(input('Is the force still wanted, 1 (Yes) or 0 (No): ')))
                    return 'False', 'False', Option, 'False', 'False', 'False', 'False'

                F_mag_oi = [F_mag1, F_mag2, F_mag3]
                d_F_oi = [d_F1, d_F2, d_F3]
                
                A = np.array([[d_F1**2, d_F1, 1],[d_F2**2, d_F2, 1],[d_F3**2, d_F3, 1]])
                B = np.array([[F_mag1],[F_mag2],[F_mag3]])
                coef = np.linalg.inv(A)@B
                a, b, c = coef[0][0], coef[1][0], coef[2][0]
                
                F_mag = (1/3)*a*(d_F3**3 - d_F1**3) + (1/2)*b*(d_F3**2 - d_F1**2) + c*(d_F3 - d_F1)
                d_F = ((1/4)*a*(d_F3**4 - d_F1**4) + (1/3)*b*(d_F3**3 - d_F1**3) + (1/2)*c*(d_F3**2 - d_F1**2))/F_mag

                fxf = (1/3)*a*X**3 + (1/2)*b*X**2 + c*X + P
                mxf = (1/4)*(1/3)*a*X**4 + (1/3)*(1/2)*b*X**3 + (1/2)*c*X**2 + P*X + Q

        else:
            print('The character entered is invalid, only P (Point Load) or D (Distributed Load) is accepted.')
            Option = bool(int(input('Is the force still wanted, 1 (Yes) or 0 (No): ')))
            return 'False', 'False', Option, 'False', 'False', 'False', 'False'

        Option = bool(int(input('Is there more Force to be added, 1 (Yes) or 0 (No): ')))
        return F_mag, d_F, Option, F_mag_oi, d_F_oi, fxf, mxf
    
    def inputMoment(self, L):
        M_mag = float(input('Enter the magnitude of the moment, with clockwise being negative (N*m): '))
        d_M = float(input('Enter the location of the moment, with left fixed end being 0 (m): '))
        if d_M > L or d_M < 0:
            print('One or more end(s) of the distributed load is located outside of the beam range.')
            Option = bool(int(input('Is the moment still wanted, 1 (Yes) or 0 (No): ')))
            return 'False', 'False', Option
        Option = bool(int(input('Is there more moment to be added, 1 (Yes) or 0 (No): ')))

        return M_mag, [d_M], Option
    
    def forceandMomentBalance(self, dS, F, dF, M, dM):
        import sympy as sp
        import numpy as np

        X = sp.symbols('X')
        Q = sp.symbols('Q')

        unknown_dict = []
        if len(dS) == 0:
            raise ValueError('The beam cannot have zero support (i.e. the beam is in motion).')
        if len(dS) != 1:
            for i in range(len(dS)):
                unknown_dict.append(sp.symbols(str(chr(i+97))))
        
        # Force Balance
        F_Balance = np.sum(F)
        for i in range(len(dS)):
            F_Balance += unknown_dict[i]
        
        # Moment Balance
        M_Balance = [np.sum(M) for i in range(len(dS))]
        for i in range(len(dS)):
            cc = 0
            while cc != len(dS):
                M_Balance[i] += unknown_dict[cc]*(dS[cc] - dS[i])
                cc += 1
            for j in range(len(F)):
                M_Balance[i] += F[j]*(dF[j] - dS[i])

        F_S = sp.solve_poly_system(M_Balance, unknown_dict)
        FS = [[i] for i in F_S[0]]
        FxS = [i for i in F_S[0]]
        MxS = [i[0]*X + Q for i in FS]
        
        return F_Balance, M_Balance, FS, FxS, MxS
    
    def drawDiagram(self, FS, dS, F_oi, dF_oi, M, dM, FxF, FxS, MxF, MxS):
        import matplotlib.pyplot as plt
        import sympy as sp
        import numpy as np

        X = sp.symbols('X')
        P = sp.symbols('P')
        Q = sp.symbols('Q')
        

        d1 = dF_oi + dS # original force and support locations only for SFD
        Fx1 = FxF + FxS # all force expressions and support expressions
        F1 = F_oi + FS # original force and force of supports in magnitude
        Mx1 = MxF + MxS # phase 1 of moment expressions, only include force and support effects, moment in magnitude is missing

        d2, Mx2 = [], []

        zip1 = sorted(zip(d1,Fx1,F1,Mx1))
        d1, Fx1, F1, Mx1 = [], [], [], []
        for i in zip1:
            d1.append(i[0])
            Fx1.append(i[1])
            F1.append(i[2])
            Mx1.append(i[3])

        points_SFD = [(0,0)]
        x, y = 0, 0
        i = 0
        removed_idx = []
        while i != len(d1):
            flag = 0
            if len(d1[i]) == 1:
                x, y = d1[i][0], y
                points_SFD.append((x,y))
                y += F1[i][0]
                points_SFD.append((x,y))
                flag = 1
            elif len(d1[i]) > 1:
                dl, dr = min(d1[i]), max(d1[i])
                d_in_between = [d1[j][0] for j in range(len(d1)) if len(d1[j]) == 1 and dl < d1[j][0] and dr > d1[j][0]]
                F_in_between = [F1[j][0] for j in range(len(d1)) if len(d1[j]) == 1 and dl < d1[j][0] and dr > d1[j][0]]
                if d_in_between != []:
                    for j in range(len(d_in_between) + 1):
                        if j == 0:
                            # dl to in_between[j][0]
                            f = sp.lambdify(X, Fx1[i], 'numpy')
                            Cc = sp.solve(f(points_SFD[-1][0]) - points_SFD[-1][1], P)
                            Fx1[i] = Fx1[i] - P + Cc[0]
                            Mx1[i] = Mx1[i] - P*X + Cc[0]*X
                            Mx2.append(Mx1[i])
                            d2.append([dl, d_in_between[j]])
                            f = sp.lambdify(X, Fx1[i], 'numpy')
                            t = np.linspace(dl, d_in_between[j], 15)
                            for k in t:
                                points_SFD.append((k,f(k)))
                            x, y = points_SFD[-1][0], points_SFD[-1][1] + F_in_between[j]
                            points_SFD.append((x,y))
                            Fx1[i] = Fx1[i] + F_in_between[j]
                            Mx1[i] = Mx1[i] + F_in_between[j]*X
                            f = sp.lambdify(X, Fx1[i], 'numpy')
                        elif j == len(d_in_between):
                            # in_between[j-1][0] to dr
                            Mx2.append(Mx1[i])
                            d2.append([d_in_between[j-1], dr])
                            t = np.linspace(d_in_between[j-1], dr, 15)
                            for k in t:
                                points_SFD.append((k,f(k)))
                            x, y = points_SFD[-1][0], points_SFD[-1][1]
                            points_SFD.append((x,y))
                        else:
                            # in_between[j-1][0] to in_between[j][0]
                            Mx2.append(Mx1[i])
                            d2.append([d_in_between[j-1], d_in_between[j]])
                            t = np.linspace(d_in_between[j-1], d_in_between[j], 15)
                            for k in t:
                                points_SFD.append((k,f(k)))
                            x, y = points_SFD[-1][0], points_SFD[-1][1] + F_in_between[j]
                            points_SFD.append((x,y))
                            Fx1[i] = Fx1[i] + F_in_between[j]
                            Mx1[i] = Mx1[i] + F_in_between[j]*X
                            f = sp.lambdify(X, Fx1[i], 'numpy')
                    
                    removed_idx.append(i)
                    i += len(d_in_between)
                    flag = 2
                elif d_in_between == []:
                    f = sp.lambdify(X, Fx1[i], 'numpy')
                    Cc = sp.solve(f(points_SFD[-1][0]) - points_SFD[-1][1], P)
                    Fx1[i] = Fx1[i] - P + Cc[0]
                    Mx1[i] = Mx1[i] - P*X + Cc[0]*X
                    f = sp.lambdify(X, Fx1[i], 'numpy')
                    t = np.linspace(dl, dr, 15)
                    for j in t:
                        points_SFD.append((j,f(j)))
                    x, y = points_SFD[-1][0], points_SFD[-1][1]
                    flag = 1
            i += 1

            if i != len(d1):
                points_SFD.append((d1[i][0],points_SFD[-1][1]))
                if flag == 1:
                    d2.append([d1[i-1][-1],d1[i][0]])
                    Mx2.append(points_SFD[-1][1]*X + Q)
                elif flag == 2:
                    d2.append([d1[i-1-len(d_in_between)][-1],d1[i][0]])
                    Mx2.append(points_SFD[-1][1]*X + Q)

        removed_idx = sorted(removed_idx, reverse = True)
        for i in removed_idx:
            del d1[i]
            del Mx1[i]

        d2 = d2 + d1 + dM
        Mx2 = Mx2 + Mx1 + [0 for i in range(len(dM))]
        M2 = [0 for i in range(len(d2) - len(dM))] + M

        zip2 = sorted(zip(d2,M2,Mx2))
        d2, Mx2, M2 = [], [], []
        for i in zip2:
            d2.append(i[0])
            M2.append(i[1])
            Mx2.append(i[2])

        # Bending Moment Diagram
        points_BMD = [(0,0)]
        x, y = 0, 0
        i = 0
        while i != len(d2):
            flag = 0
            if len(d2[i]) == 1:
                x, y = d2[i][0], y
                points_BMD.append((x,y))
                y += M2[i]
                points_BMD.append((x,y))
                flag = 1
            elif len(d2[i]) > 1:
                dl, dr = min(d2[i]), max(d2[i])
                d_in_between = [d2[j][0] for j in range(len(d2)) if len(d2[j]) == 1 and dl < d2[j][0] and dr > d2[j][0]]
                M_in_between = [M2[j] for j in range(len(d2)) if len(d2[j]) == 1 and dl < d2[j][0] and dr > d2[j][0]]
                if d_in_between != []:
                    for j in range(len(d_in_between) + 1):
                        if j == 0:
                            # dl to in_between[j][0]
                            f = sp.lambdify(X, Mx2[i], 'numpy')
                            Cc = sp.solve(f(points_BMD[-1][0]) - points_BMD[-1][1], Q)
                            if Cc == []:
                                Cc = [0]
                            Mx2[i] = Mx2[i] - Q + Cc[0]
                            f = sp.lambdify(X, Mx2[i], 'numpy')
                            t = np.linspace(dl, d_in_between[j], 15)
                            for k in t:
                                points_BMD.append((k,f(k)))
                            x, y = points_BMD[-1][0], points_BMD[-1][1] + M_in_between[j]
                            points_BMD.append((x,y))
                            Mx2[i] = Mx2[i] + M_in_between[j]
                            f = sp.lambdify(X, Mx2[i], 'numpy')
                        elif j == len(d_in_between):
                            # in_between[j-1][0] to dr
                            t = np.linspace(d_in_between[j-1], dr, 15)
                            for k in t:
                                points_BMD.append((k,f(k)))
                            x, y = points_BMD[-1][0], points_BMD[-1][1]
                            points_BMD.append((x,y))
                        else:
                            # in_between[j-1][0] to in_between[j][0]
                            t = np.linspace(d_in_between[j-1], d_in_between[j], 15)
                            for k in t:
                                points_BMD.append((k,f(k)))
                            x, y = points_BMD[-1][0], points_BMD[-1][1] + M_in_between[j]
                            points_BMD.append((x,y))
                            Mx2[i] = Mx2[i] + M_in_between[j]
                            f = sp.lambdify(X, Mx2[i], 'numpy')

                    i += len(d_in_between)
                    flag = 2
                elif d_in_between == []:
                    f = sp.lambdify(X, Mx2[i], 'numpy')
                    Cc = sp.solve(f(points_BMD[-1][0]) - points_BMD[-1][1], Q)
                    Mx2[i] = Mx2[i] - Q + Cc[0]
                    f = sp.lambdify(X, Mx2[i], 'numpy')
                    t = np.linspace(dl, dr, 15)
                    for j in t:
                        points_BMD.append((j,f(j)))
                    x, y = points_BMD[-1][0], points_BMD[-1][1]
                    flag = 1
            i += 1



        points_SFD_x = [i[0] for i in points_SFD]
        points_SFD_y = [i[1] for i in points_SFD]
        # print(points_SFD)
        points_BMD_x = [i[0] for i in points_BMD]
        points_BMD_y = [i[1] for i in points_BMD]
        # print(points_BMD)


        plot, (ax1, ax2) = plt.subplots(2,1)
        ax1.plot(points_SFD_x, points_SFD_y)
        ax1.grid(True, which = 'both')

        ax1.axhline(y = 0, color = 'k')
        ax1.axvline(x = 0, color = 'k')

        ax2.plot(points_BMD_x, points_BMD_y)
        ax2.grid(True, which = 'both')

        ax2.axhline(y = 0, color = 'k')
        ax2.axvline(x = 0, color = 'k')

        plt.show()

    def main(self):
        fg1, fg2 = True, True
        L = self.inputPara()
        dS = self.inputSupport(L)
        dS = sorted(dS)
        print('With left end being 0 m, the location of two supports are:', dS)
        F, dF = [],[]
        M, dM = [],[]
        F_oi, dF_oi = [],[]
        FxF, MxF = [], []

        while fg1 == True:
            F_mag, d_F, fg1, F_mag_oi, d_F_oi, fxf, mxf = self.inputForce(L)
            if F_mag != 'False' and d_F != 'False' and F_mag_oi != 'False' and d_F_oi != 'False' and fxf != 'False' and mxf != 'False':
                F.append(F_mag)
                dF.append(d_F)
                F_oi.append(F_mag_oi)
                dF_oi.append(d_F_oi)
                FxF.append(fxf)
                MxF.append(mxf)
        fg2 = bool(int(input('Is there a point moment present, 1 (Yes) or 0 (No): ')))
        while fg2 == True:
            M_mag, d_M, fg2 = self.inputMoment(L)
            if M_mag != 'False' and d_M != 'False':
                M.append(M_mag)
                dM.append(d_M)

        print('The Force(s) are ', F, 'at the location of ', dF)
        if len(M) != 0:
            print('The Moment(s) are ', M, 'at the location of ', dM)
        else:
            print('There is no Moment present in on the given beam.')

        F_Balance, M_Balance, FS, FxS, MxS = self.forceandMomentBalance(dS, F, dF, M, dM)
        dS = [[i] for i in dS]
        M = [-i for i in M]
        print('Force Equalibrium Equation is ', F_Balance)
        print('Moment Equalibrium Equations regarding two supports are ', M_Balance)

        # print('FS = ', FS)
        # print('dS = ', dS)
        # print('F_oi = ', F_oi)
        # print('dF_oi = ', dF_oi)
        # print('M = ', M)
        # print('dM = ', dM)
        # print('FxF = ', FxF)
        # print('FxS = ', FxS)
        # print('MxF = ', MxF)
        # print('MxS = ', MxS)

        self.drawDiagram(FS, dS, F_oi, dF_oi, M, dM, FxF, FxS, MxF, MxS)

calcSFDandBMD = SFDandBMD()
calcSFDandBMD.main()