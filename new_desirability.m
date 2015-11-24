function Phi=new_desirability(h,theta,t,K,alpha)

Beta = K*(t+alpha);

phi  = bsxfun(@times,exp(-Beta*theta),h);
Phi  = (t+alpha)*trapz(theta,phi);

%plot(t,h,t,phi)