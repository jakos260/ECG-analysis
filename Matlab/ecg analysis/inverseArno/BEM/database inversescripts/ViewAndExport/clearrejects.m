%clearrejects
global bs_Points bs_PointType bs_rejected
bs_Points(bs_rejected==1,:)=0;
bs_PointType(bs_rejected==1,:)=0;