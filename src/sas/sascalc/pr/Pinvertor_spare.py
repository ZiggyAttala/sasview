    def get_matrix(self, nfunc, nr, a_obj, b_obj):
        """
        Returns A matrix and b vector for least square problem.
        @param nfunc: number of base functions
        @param nr: number of r-points used when evaluating reg term.
        @param a: A array to fill
        @param b: b vector to fill
        @return: 0
        """
        #replace assert(n_b>=nfunc) and assert(n_a>=nfunc*(nr+self.npoints))

        if(b_obj.shape[0] < nfunc):
            #Vector too small
            #abort program execution, replace C assert()
            #may need more error messages here etc.
            print("b too small")
            return None
        if(a_obj.shape[0] < nfunc*(nr + self.npoints)):
            print("a too small")
            return None

        a = a_obj
        b = b_obj

        sqrt_alpha = np.sqrt(self.alpha)
        pi = np.arccos(-1.0)
        offset = (1, 0)[self.est_bck == 1]

        #instead of checking for 0 in err in for loop, check all
        #for 0 before
        def check_for_zero(x):
            for i, ni in enumerate(x):
                if(ni == 0):
                    return True
            return False

        if(check_for_zero(self.err)):
            logger.error("Pinvertor.get_matrix: Some I(Q) points have no error.")
            return None

        for j in range(nfunc):
            for i in range(self.npoints):
                index = i * nfunc + j
                npts = 21
                if(self.accept_q(self.x[i])):

                    if(self.est_bck == 1 and j == 0):
                        a[index] = 1.0/self.err[i]

                    else:
                        if(self.slit_width > 0 or self.slit_height > 0):
                            a[index] = py_invertor.ortho_transformed_smeared_qvec_njit(self.x[i], self.d_max, j + offset,
                                                                             self.slit_height, self.slit_width, npts)/self.err[i]
                        else:
                            a[index] = py_invertor.ortho_transformed_qvec_njit(self.x[i], self.d_max, j + offset)/self.err[i]

            for i_r in range(nr):
                index = (i_r + self.npoints) * nfunc + j
                if(self.est_bck == 1 and j == 0):
                    a[index] = 0.0
                else:
                    r = self.d_max / nr * i_r
                    tmp = pi * (j + offset) / self.d_max
                    t1 = sqrt_alpha * 1.0/nr * self.d_max * 2.0

                    t2 = (2.0 * pi * (j + offset)/self.d_max * np.cos(pi * (j + offset)*r/self.d_max)
                    + tmp * tmp * r * np.sin(pi * (j + offset)*r/self.d_max))

                    a[index] =  t1 * t2

        for i in range(self.npoints):
            if(self.accept_q(self.x[i])):
                b[i] = self.y[i] / self.err[i]

        return 0