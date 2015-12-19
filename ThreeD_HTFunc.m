function [ Mapn ] = ThreeDFiniteDiff( Map,RCP,kmat,ts,ds,Ta,dz,hs)
%ThreeDFiniteDiff
%This function takes in parameters and returns a heatmap of an area at a
%certain time.  Evaluate in a loop to step through times.

[rw,cl,bx]=size(Map);
Mapn = zeros(rw,cl,bx);
%Set up functions.

mxyz = ts/(RCP*ds^2*dz);


%Use finite difference equations to get the new heat map.

%There's eight corners.  A corner is all three extremities.
for a = 1:bx
    for k = 1:rw
        for n = 1 : cl
        
            
            if k == 1
                if n == 1
                    %Two Corners here and fixed.
                    %Bottom left edge
                    if a == 1 %(0,0,0)
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n+1,a)+Map(k+1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))))+ Map(k,n,a);
                    elseif a == bx %(0,0,H)
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+hs*ds^2/4*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n+1,a)+Map(k+1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))))+ Map(k,n,a);
                    else %(0,y,H)
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat/4*(2*dz*(Map(k,n+1,a)+Map(k+1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)+Map(k,n,a+1)-2*Map(k,n,a))))+ Map(k,n,a);
                    end
                    
                elseif n==cl %Two more corners here.
                    %Bottom right edge
                    if a == 1
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n-1,a)+Map(k+1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))))+ Map(k,n,a);
                    elseif a == bx
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+hs*ds^2/4*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n-1,a)+Map(k+1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))))+ Map(k,n,a);
                    else
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat/4*(2*dz*(Map(k,n-1,a)+Map(k+1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)+Map(k,n,a+1)-2*Map(k,n,a))))+ Map(k,n,a);
                    end
                    
                else
                    %Bottom face
                    if a == 1
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a)) + kmat/4*(dz*(Map(k,n+1,a)+Map(k,n-1,a)-2*Map(k,n,a))+2*ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))+2*dz*(Map(k+1,n,a)-Map(k,n,a))))+ Map(k,n,a);
                    elseif a == bx
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+hs*ds^2/2*(Ta-Map(k,n,a)) + kmat/4*(2*dz*(Map(k+1,n,a)-Map(k,n,a))+2*ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))+dz*ds*(Map(k,n-1,a)+Map(k,n+1,a)-2*Map(k,n,a))))+ Map(k,n,a);
                    else %(0,y,z)
                        Mapn(k,n,a) = 2*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat*(ds^2/(2*dz)*(Map(k,n,a+1)+Map(k,n,a-1)-2*Map(k,n,a)) + dz/2*(Map(k,n+1,a)+Map(k,n-1,a)-2*Map(k,n,a))+dz*(Map(k+1,n,a)-Map(k,n,a)))) + Map(k,n,a);
                    end
                end
                
            elseif k == rw
                
                if n == 1
                    
                    %Two corners here and Fixed
                    %Top left edge
                    if a == 1
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n+1,a)+Map(k-1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))))+Map(k,n,a);
                    elseif a == bx
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+hs*ds^2/4*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n+1,a)+Map(k-1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))))+Map(k,n,a);
                    else
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat/4*(2*dz*(Map(k,n+1,a)+Map(k-1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)+Map(k,n,a+1)-2*Map(k,n,a))))+ Map(k,n,a);
                    end
                    %Other two corners here
                    %Top right edge
                elseif n == cl
                    if a == 1
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n-1,a)+Map(k-1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))))+Map(k,n,a);
                    elseif a == bx
                        Mapn(k,n,a) = 8*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a))+hs*ds^2/4*(Ta-Map(k,n,a))+kmat/4*(dz*(Map(k,n-1,a)+Map(k-1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))))+Map(k,n,a);
                    else
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat/4*(2*dz*(Map(k,n-1,a)+Map(k-1,n,a)-2*Map(k,n,a))+ds^2/dz*(Map(k,n,a-1)+Map(k,n,a+1)-2*Map(k,n,a))))+ Map(k,n,a);
                    end
                    
                    %Top face
                else
                    if a == 1
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a)) + kmat/4*(dz*(Map(k,n+1,a)+Map(k,n-1,a)-2*Map(k,n,a))+2*ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))+2*dz*(Map(k-1,n,a)-Map(k,n,a))))+ Map(k,n,a);
                    elseif a == bx
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a)) + hs*ds^2/2*(Ta-Map(k,n,a)) + kmat/4*(2*dz*(Map(k-1,n,a)-Map(k,n,a))+2*ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))+dz*ds*(Map(k,n-1,a)+Map(k,n+1,a)-2*Map(k,n,a))))+ Map(k,n,a);
                    else %(L,y,z)
                        Mapn(k,n,a) = 2*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat*(ds^2/(2*dz)*(Map(k,n,a+1)+Map(k,n,a-1)-2*Map(k,n,a)) + dz/2*(Map(k,n+1,a)+Map(k,n-1,a)-2*Map(k,n,a))+dz*(Map(k-1,n,a)-Map(k,n,a)))) + Map(k,n,a);
                    end
                    
                end
                
            else
                
                if n == 1
                    
                    %Left Face
                    if a == 1
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a)) + kmat/4*(dz*(Map(k+1,n,a)+Map(k-1,n,a)-2*Map(k,n,a))+2*ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))+2*dz*(Map(k,n+1,a)-Map(k,n,a))))+ Map(k,n,a);
                    elseif a == bx
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a)) + hs*ds^2/2*(Ta-Map(k,n,a)) + kmat/4*(dz*(Map(k+1,n,a)+Map(k-1,n,a)-2*Map(k,n,a))+2*ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))+2*dz*(Map(k,n+1,a)-Map(k,n,a))))+ Map(k,n,a);
                    else
                        Mapn(k,n,a) = 2*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat*(ds^2/(2*dz)*(Map(k,n,a+1)+Map(k,n,a-1)-2*Map(k,n,a))+dz/2*(Map(k+1,n,a)+Map(k-1,n,a)-2*Map(k,n,a))+dz*(Map(k,n+1,a)-Map(k,n,a))))+Map(k,n,a);
                    end
                    
                elseif n == cl
                    
                    %Right Face
                    if a == 1
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a)) + kmat/4*(dz*(Map(k+1,n,a)+Map(k-1,n,a)-2*Map(k,n,a))+2*ds^2/dz*(Map(k,n,a+1)-Map(k,n,a))+2*dz*(Map(k,n-1,a)-Map(k,n,a))))+ Map(k,n,a);
                    elseif a == bx
                        Mapn(k,n,a) = 4*mxyz*(hs*ds*dz/2*(Ta-Map(k,n,a)) + hs*ds^2/2*(Ta-Map(k,n,a)) + kmat/4*(dz*(Map(k+1,n,a)+Map(k-1,n,a)-2*Map(k,n,a))+2*ds^2/dz*(Map(k,n,a-1)-Map(k,n,a))+2*dz*(Map(k,n-1,a)-Map(k,n,a))))+ Map(k,n,a);
                    else
                        Mapn(k,n,a) = 2*mxyz*(hs*ds*dz*(Ta-Map(k,n,a)) + kmat*(ds^2/(2*dz)*(Map(k,n,a+1)+Map(k,n,a-1)-2*Map(k,n,a))+dz/2*(Map(k+1,n,a)+Map(k-1,n,a)-2*Map(k,n,a))+dz*(Map(k,n-1,a)-Map(k,n,a))))+Map(k,n,a);
                    end
                    
                else
                    
                    if a == 1
                        %Underside (OK)
                        Mapn(k,n,a) = 2*mxyz*kmat*(dz/2*(Map(k,n+1,a) + Map(k,n-1,a) + Map(k+1,n,a) + Map(k-1,n,a)-4*Map(k,n,a)) + ds^2/dz *(Map(k,n,a+1)-Map(k,n,a)))+ Map(k,n,a);
                    elseif a == bx
                        %Topside
                        Mapn(k,n,a) = 2*mxyz*(hs*ds^2*(Ta-Map(k,n,a))+kmat*(dz*(Map(k,n+1,a) + Map(k,n-1,a) + Map(k+1,n,a) + Map(k-1,n,a)-4*Map(k,n,a)) + ds^2/dz *(Map(k,n,a-1)-Map(k,n,a))))+ Map(k,n,a);
                    else
                        %Center (OK)
                        Mapn(k,n,a) = mxyz*kmat*(dz*(Map(k,n+1,a) + Map(k,n-1,a) + Map(k+1,n,a) + Map(k-1,n,a)-4*Map(k,n,a)) + ds^2/dz *(Map(k,n,a+1)-2*Map(k,n,a)+Map(k,n,a-1)))+ Map(k,n,a);
                    end
                end
            end
        end
    end
end


end

