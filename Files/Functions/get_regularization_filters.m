function [ Hreg ] = get_regularization_filters( w, fc )


   % fc = interp1((1:4)',[90,680,1650,2600]'/10,(1:Nmax),'linear','extrap');
fc(fc>20e3)=20e3;

f = w / (2*pi);
nn = (2:length(fc)+1);

Hl = 1./(1+(f./fc).^nn);
Hu = (f./fc).^nn./(1+(f./fc).^nn);
H = [Hl,ones(size(Hl,1),1)].*[ones(size(Hl,1),1),Hu];
Hs = sum(H,2);
H = H./Hs;
Hreg =fliplr( cumsum(fliplr(H),2) );
Hreg(isnan(Hreg)) = 0;


end

