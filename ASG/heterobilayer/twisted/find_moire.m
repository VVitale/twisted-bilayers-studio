function [m_s,n_s,p_s,theta_s,nmoire] = find_moire(theta,dt,r_t,max_int,size_t,comb,threshold,p,r_p,m_s,n_s,theta_s,p_s)
         deg2rad=pi/180;
         rad2deg = 1.0/deg2rad;
         
         dis_m = -1;
         dis_n = -1;
         xs = -99999;
         
         iter = 0;
         nmoire = 0;
         
         % Allocate arrays
         x_m = zeros(2,size_t);
         x_n = zeros(2,size_t);
         x_m(:,:) = -88888;
         x_n(:,:) = -99999;
         lfound = false(size_t,1);

         % Start search
         for th = theta - r_t : dt : theta + r_t
             iter = iter + 1;
             th = th * deg2rad;
             ct = cos(th);
             st = sin(th);
             for im = 1 : max_int
                 if (lfound(iter))
                     break;
                 end
                 dis_m = (ct + 1/(2*im)*(st/sqrt(3)+ct))^2 - 1.0 - 1.0/im;
                 if (dis_m >= 0)
                     x_m(1,iter) = ct + 1/(2*im)*(st/sqrt(3) + ct) + sqrt(dis_m);
                     x_m(2,iter) = ct + 1/(2*im)*(st/sqrt(3) + ct) - sqrt(dis_m);
                 end
                 
                 if (comb == -1)
                     start_in = im + 1;
                 else
                     start_in = 1;
                 end
                 for in = start_in : max_int
                     if (lfound(iter) == true)
                         break;
                     end
                     xs = -99999;
                     err = -99999;
                     ierr = 0;
                     dis_n = (ct+st/(sqrt(3)*in))^2-1;
                     if (dis_n >= 0)
                         x_n(1,iter) = ct + 1/(in)*st/sqrt(3) + sqrt(dis_n);
                         x_n(2,iter) = ct + 1/(in)*st/sqrt(3) - sqrt(dis_n);
                     end
                     % The assumption here is that we have only one solution
                     if (abs(x_m(1,iter) - x_n(1,iter)) < threshold)
                         err = abs(x_m(1,iter) - x_n(1,iter));
                         xs = x_m(1,iter);
                     elseif (abs(x_m(1,iter) - x_n(2,iter)) < threshold)
                         err = abs(x_m(1,iter) - x_n(2,iter));
                         xs = x_m(1,iter);                              
                     elseif (abs(x_m(2,iter) - x_n(1,iter)) < threshold)
                         err = abs(x_m(2,iter) - x_n(1,iter));
                         xs = x_m(2,iter);
                     elseif (abs(x_m(2,iter) - x_n(2,iter)) < threshold)
                         err = abs(x_m(2,iter) - x_n(2,iter));
                         xs = x_m(2,iter);
                     end
                     if (err < threshold)
                         if (abs(xs - p) <= r_p)
                             % Check that we have not found it yet
                             for i = 1 : iter - 1;
                                 % if found
                                 if (im == m_s(i)) && (in == n_s(i))
                                     lfound(iter) = true;
                                     loc_err = -99999;
                                     msg = ['Found another structure with m = ', ...
                                         num2str(im), '; n = ', num2str(in), ...
                                         ' (theta = ', num2str(th*rad2deg),  ...
                                         num2str(in), '; p = ', num2str(xs), ... 
                                         '; err = ', num2str(err),')', '...'];
                                     disp(msg)
                                     % Check if the error is smaller and swap if it
                                     % is
                                     loc_err = min([abs(x_m(1,i)-x_n(1,i)),abs(x_m(1,i)-x_n(2,i)),...
                                             abs(x_m(2,i)-x_n(1,i)),abs(x_m(2,i)-x_n(2,i))]);
                                     if (loc_err > err)
                                         m_s(i) = im;
                                         n_s(i) = in;
                                         theta_s(i) = th;
                                         p_s(i) = xs;
                                         x_m(:,i) = x_m(:,iter);
                                         x_n(:,i) = x_n(:,iter);
                                         msg = 'Accepted';
                                         disp(msg)
                                     else
                                         msg = 'Rejected';
                                         disp(msg)
                                     end
                                 end
                             end
                             if (~lfound(iter))
                                     lfound(iter) = true;
                                     nmoire = nmoire + 1;
                                     m_s(iter) = im;
                                     n_s(iter) = in;
                                     theta_s(iter) = th;
                                     p_s(iter) = xs;
                                     msg = ['Found one moire supercell with theta = ', ...
                                         num2str(th*rad2deg),'; m = ',num2str(im), ...
                                         '; n = ', num2str(in), '; p = ', num2str(xs),...
                                         '; err = ', num2str(err)];
                                     disp(msg)
                                     break;
                             end
                         end
                     end
                 end
             end
         end
end
