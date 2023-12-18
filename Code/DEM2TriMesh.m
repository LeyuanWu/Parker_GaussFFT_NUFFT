function [Vert,Faces]=DEM2TriMesh(HDEM,href,nEst,nNor,Estmin,Estmax,Normin,Normax,eswn,flagTri)
% ----------------------------------------------------------------------- %
% M-file to convert a DEM to a triangular mesh polyhedral model
% ----------------------------------------------------------------------- %
% Note: the coordinate system of east, north and upward is adopted for common DEM data
% ---- input ---- %
% HDEM: DEM data
% href: Reference height
% nEst,nNor,Estmin,Estmax,Normin,Normax: DEM Grid size and range
% eswn: flag indicating the triangular mesh direction
%       'ne': southwest-northeast triangulation
%       'nw': southeast-northwest triangulation
% flagTri: 0 The terrain is triangular, the sides and bottom are polygonal
%          1 All triangular, including the sides and bottom
%          2 Triangular only for the terrain
% ---- output ---- %
% Vert,Faces: Vertex and face lists
% ---- Note ---- %
% When generating polygons with non-fully triangular faces, it is very 
% memory consuming when the number of edges for the side faces are large.
% Using cell array instead would be much better.
%%%%%%%% DEM grid generation
Est=linspace(Estmin,Estmax,nEst);Est=Est';
Nor=linspace(Normin,Normax,nNor);Nor=Nor';
[EST,NOR]=meshgrid(Est,Nor);
szGrid=size(EST);
%%%%%%%% Number of grid cells
mEst=nEst-1;mNor=nNor-1;
%%%%%%%% Combination of all indices
[indEST,indNOR]=meshgrid(1:1:nEst,1:1:nNor);
tempEst1=indEST(1:end-1,1:end-1);
tempEst2=indEST(2:end,1:end-1);
tempEst3=indEST(2:end,2:end);
tempEst4=indEST(1:end-1,2:end);
tempNor1=indNOR(1:end-1,1:end-1);
tempNor2=indNOR(2:end,1:end-1);
tempNor3=indNOR(2:end,2:end);
tempNor4=indNOR(1:end-1,2:end);
%%%%%%%% Vertex list
VertEst=EST(:);VertNor=NOR(:);VertH=HDEM(:); 
N=nEst*nNor;
VertB1=[Estmin,Normin,href];VertB2=[Estmin,Normax,href];
VertB3=[Estmax,Normax,href];VertB4=[Estmax,Normin,href];
% ----------------------------------------------------------------------- %
% Add vertices to the four sides of the bottom plane.
%                 2  2  2  2  2  2  2  2  3
%                 1                       3
%                 1                       3
%                 1                       3
%                 1                       3
%                 1  4  4  4  4  4  4  4  4
% ----------------------------------------------------------------------- %
HREF1=[Estmin*ones(mNor,1),Nor(1:end-1),href*ones(mNor,1)]; 
HREF2=[Est(1:end-1),Normax*ones(mEst,1),href*ones(mEst,1)];
HREF3=[Estmax*ones(mNor,1),Nor(end:-1:2),href*ones(mNor,1)];
HREF4=[Est(end:-1:2),Normin*ones(mEst,1),href*ones(mEst,1)];
if(flagTri==0)
    Vert=vertcat([VertEst,VertNor,VertH],VertB1,VertB2,VertB3,VertB4);
end
if(flagTri==1)
    Vert=vertcat([VertEst,VertNor,VertH],HREF1,HREF2,HREF3,HREF4);
end
if(flagTri==2)
    Vert=vertcat([VertEst,VertNor,VertH]);
end
%%%%%%%% Face list
if(flagTri==0)
    nFace=mEst*mNor*2+4+1;
    Faces=cell(nFace,1);
    if(strcmp(eswn,'ne'))
        tempUL=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        tempDR=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)));
        Faces(1:mEst*mNor)=mat2cell(tempUL,ones(mEst*mNor,1));
        Faces(mEst*mNor+1:2*mEst*mNor)=mat2cell(tempDR,ones(mEst*mNor,1));
    end
    if(strcmp(eswn,'nw'))
        tempUL=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        tempDR=horzcat(sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        Faces(1:mEst*mNor)=mat2cell(tempUL,ones(mEst*mNor,1));
        Faces(mEst*mNor+1:2*mEst*mNor)=mat2cell(tempDR,ones(mEst*mNor,1));
    end
    %%%% left side
    lFace=sub2ind(szGrid, 1:1:nNor, ones(1,nNor));
    lFace=horzcat(lFace,nEst*nNor+2,nEst*nNor+1);
    Faces{mEst*mNor*2+1}=lFace;
    %%%% right side
    rFace=sub2ind(szGrid, nNor:-1:1, nEst*ones(1,nNor));
    rFace=horzcat(rFace,nEst*nNor+4,nEst*nNor+3);
    Faces{mEst*mNor*2+2}=rFace;
    %%%% front side
    fFace=sub2ind(szGrid, nNor*ones(1,nEst), 1:1:nEst);
    fFace=horzcat(fFace,nEst*nNor+3,nEst*nNor+2);
    Faces{mEst*mNor*2+3}=fFace;
    %%%% back side
    bFace=sub2ind(szGrid, ones(1,nEst), nEst:-1:1);
    bFace=horzcat(bFace,nEst*nNor+1,nEst*nNor+4);
    Faces{mEst*mNor*2+4}=bFace;
    %%%% bottom plane
    btm=[nEst*nNor+1,nEst*nNor+2,nEst*nNor+3,nEst*nNor+4];
    Faces{mEst*mNor*2+5}=btm;
end
if(flagTri==1)
    nFace=mEst*mNor*2+2*mNor+2*mEst+2*mNor+2*mEst+2;
    Faces=zeros(nFace,3);
    if(strcmp(eswn,'ne'))
        tempUL=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        tempDR=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)));
        Faces(1:mEst*mNor,:)=tempUL;
        Faces(mEst*mNor+1:2*mEst*mNor,:)=tempDR;
    end
    if(strcmp(eswn,'nw'))
        tempUL=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        tempDR=horzcat(sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        Faces(1:mEst*mNor,:)=tempUL;
        Faces(mEst*mNor+1:2*mEst*mNor,:)=tempDR;
    end
    %%%% left side
    for i=1:1:mNor
        Faces(mEst*mNor*2+2*(i-1)+1,:)=[N+i,sub2ind(szGrid,i,1),N+i+1];
        Faces(mEst*mNor*2+2*(i-1)+2,:)=[N+i+1,sub2ind(szGrid,i,1),sub2ind(szGrid,i+1,1)];
    end
    %%%% front side
    for i=1:1:mEst
        Faces(mEst*mNor*2+2*mNor+2*(i-1)+1,:)=[N+mNor+i,sub2ind(szGrid,nNor,i),N+mNor+i+1];
        Faces(mEst*mNor*2+2*mNor+2*(i-1)+2,:)=[N+mNor+i+1,sub2ind(szGrid,nNor,i),sub2ind(szGrid,nNor,i+1)];
    end
    %%%% right side
    for i=1:1:mNor
        Faces(mEst*mNor*2+2*mNor+2*mEst+2*(i-1)+1,:)=[N+mNor+mEst+i,sub2ind(szGrid,nNor+1-i,nEst),N+mNor+mEst+i+1];
        Faces(mEst*mNor*2+2*mNor+2*mEst+2*(i-1)+2,:)=[N+mNor+mEst+i+1,sub2ind(szGrid,nNor+1-i,nEst),sub2ind(szGrid,nNor-i,nEst)];
    end
    %%%% back side
    for i=1:1:mEst-1
        Faces(mEst*mNor*2+2*mNor+2*mEst+2*mNor+2*(i-1)+1,:)=[N+mNor+mEst+mNor+i,sub2ind(szGrid,1,nEst+1-i),N+mNor+mEst+mNor+i+1];
        Faces(mEst*mNor*2+2*mNor+2*mEst+2*mNor+2*(i-1)+2,:)=[N+mNor+mEst+mNor+i+1,sub2ind(szGrid,1,nEst+1-i),sub2ind(szGrid,1,nEst-i)];
    end
    Faces(mEst*mNor*2+2*mNor+2*mEst+2*mNor+2*(mEst-1)+1,:)=[N+mNor+mEst+mNor+mEst,sub2ind(szGrid,1,2),N+1];
    Faces(mEst*mNor*2+2*mNor+2*mEst+2*mNor+2*(mEst-1)+2,:)=[N+1,sub2ind(szGrid,1,2),sub2ind(szGrid,1,1)];
    %%%% bottom plane
    Faces(mEst*mNor*2+2*mNor+2*mEst+2*mNor+2*mEst+1,:)=[N+1,N+mNor+1,N+mNor+mEst+1];
    Faces(mEst*mNor*2+2*mNor+2*mEst+2*mNor+2*mEst+2,:)=[N+1,N+mNor+mEst+1,N+mNor+mEst+mNor+1];
end
if(flagTri==2)
    nFace=mEst*mNor*2;
    Faces=zeros(nFace,3);
    if(strcmp(eswn,'ne'))
        tempUL=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        tempDR=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)));
        Faces(1:mEst*mNor,:)=tempUL;
        Faces(mEst*mNor+1:2*mEst*mNor,:)=tempDR;
    end
    if(strcmp(eswn,'nw'))
        tempUL=horzcat(sub2ind(szGrid,tempNor1(:),tempEst1(:)),sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        tempDR=horzcat(sub2ind(szGrid,tempNor4(:),tempEst4(:)),sub2ind(szGrid,tempNor3(:),tempEst3(:)),sub2ind(szGrid,tempNor2(:),tempEst2(:)));
        Faces(1:mEst*mNor,:)=tempUL;
        Faces(mEst*mNor+1:2*mEst*mNor,:)=tempDR;
    end
end
end

