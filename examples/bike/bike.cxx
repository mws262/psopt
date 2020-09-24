#include "psopt.h"

const adouble m = 32.1;
const adouble g = 9.81;
const adouble Ixx = 0;
const adouble Iyy = 0;
const adouble Izz = 0;
const adouble Ixz = 0;
const adouble wheelbase = 1.05;
const adouble com_to_rear = 0.546;
const adouble com_height = 0.58;
const adouble rwheel_radius = 0.272;
const adouble fwheel_radius = 0.272;
const adouble front_wheel_kT = 1.2975;
const adouble rear_wheel_kT = 0.956;


adouble psid_eqn(adouble& phi, adouble& V, adouble& delta) {
    return -(V*sin(delta))/(wheelbase*cos(delta)*cos(phi));
}

adouble Vd_eqn(adouble& phi, adouble& V, adouble& phid, adouble& delta, adouble& deltad, adouble& torque_front, adouble& torque_rear) {
  adouble t2 = cos(delta);
  adouble t3 = cos(phi);
  adouble t4 = sin(delta);
  adouble t5 = sin(phi);
  adouble t6 = Ixz*Ixz;
  adouble t7 = V*V;
  adouble t8 = com_height*com_height;
  adouble t9 = com_height*com_height*com_height;
  adouble t11 = com_to_rear*com_to_rear;
  adouble t12 = m*m;
  adouble t13 = phid*phid;
  adouble t14 = wheelbase*wheelbase;
  adouble t15 = wheelbase*wheelbase*wheelbase;
  adouble t10 = t8*t8;
  adouble t16 = t2*t2;
  adouble t17 = t2*t2*t2;
  adouble t19 = t3*t3;
  adouble t20 = t3*t3*t3;
  adouble t22 = t4*t4;
  adouble t23 = t4*t4*t4;
  adouble t24 = t5*t5;
  adouble t25 = t5*t5*t5;
  adouble t18 = t16*t16;
  adouble t21 = t19*t19;
  adouble t26 = -t16;
  adouble t27 = t16*t19;
  adouble t28 = t26+t27+1.0;
  adouble t29 = 1.0/sqrt(t28);
 return (Ixx*fwheel_radius*t15*t17*t20*torque_rear+Ixx*rwheel_radius*t15*t27*t29*torque_front-Ixx*rwheel_radius*t15*t18*t19*t29*torque_front+Ixx*rwheel_radius*t15*t18*t21*t29*torque_front+fwheel_radius*m*t8*t15*t17*t20*torque_rear+m*rwheel_radius*t8*t15*t27*t29*torque_front-Ixz*Iyy*fwheel_radius*rwheel_radius*t5*t7*t19*t23+Ixz*Izz*fwheel_radius*rwheel_radius*t5*t7*t19*t23+V*deltad*fwheel_radius*rwheel_radius*t6*t20*t23*wheelbase-m*rwheel_radius*t8*t15*t18*t19*t29*torque_front+m*rwheel_radius*t8*t15*t18*t21*t29*torque_front+Ixx*Ixz*fwheel_radius*rwheel_radius*t4*t5*t13*t14*t27-Ixz*fwheel_radius*m*rwheel_radius*t5*t7*t8*t19*t23+V*deltad*fwheel_radius*rwheel_radius*t4*t6*t16*t20*wheelbase-com_to_rear*fwheel_radius*rwheel_radius*t5*t7*t9*t12*t19*t23-Ixx*Iyy*V*deltad*fwheel_radius*rwheel_radius*t3*t23*wheelbase+Ixx*Iyy*V*deltad*fwheel_radius*rwheel_radius*t20*t23*wheelbase-Ixx*Izz*V*deltad*fwheel_radius*rwheel_radius*t20*t23*wheelbase+Ixx*Iyy*V*deltad*fwheel_radius*rwheel_radius*t3*t4*t26*wheelbase+Ixx*Iyy*V*deltad*fwheel_radius*rwheel_radius*t4*t16*t20*wheelbase+Ixx*Izz*V*deltad*fwheel_radius*rwheel_radius*t4*t20*t26*wheelbase-Ixx*Iyy*V*fwheel_radius*phid*rwheel_radius*t2*t5*t22*wheelbase-Ixx*V*deltad*fwheel_radius*m*rwheel_radius*t3*t11*t23*wheelbase-Iyy*V*deltad*fwheel_radius*m*rwheel_radius*t3*t8*t23*wheelbase+Iyy*V*deltad*fwheel_radius*m*rwheel_radius*t8*t20*t23*wheelbase-Izz*V*deltad*fwheel_radius*m*rwheel_radius*t8*t20*t23*wheelbase-Iyy*com_height*com_to_rear*fwheel_radius*m*rwheel_radius*t5*t7*t19*t23+Izz*com_height*com_to_rear*fwheel_radius*m*rwheel_radius*t5*t7*t19*t23-Ixx*V*fwheel_radius*m*phid*rwheel_radius*t2*t5*t11*wheelbase+Ixx*V*fwheel_radius*m*phid*rwheel_radius*t5*t11*t17*wheelbase+Ixz*com_height*fwheel_radius*m*rwheel_radius*t2*t7*t20*t22*wheelbase+Ixz*fwheel_radius*m*rwheel_radius*t4*t5*t8*t13*t14*t27+V*deltad*fwheel_radius*rwheel_radius*t5*t9*t12*t14*t17*t19+V*deltad*fwheel_radius*rwheel_radius*t4*t8*t11*t12*t20*wheelbase-V*deltad*fwheel_radius*rwheel_radius*t3*t8*t11*t12*t23*wheelbase-V*deltad*fwheel_radius*rwheel_radius*t3*t10*t12*t23*t24*wheelbase+V*fwheel_radius*phid*rwheel_radius*t4*t9*t12*t14*t16*t20*2.0-V*fwheel_radius*phid*rwheel_radius*t2*t5*t8*t11*t12*wheelbase+V*fwheel_radius*phid*rwheel_radius*t5*t8*t11*t12*t17*wheelbase+V*fwheel_radius*phid*rwheel_radius*t2*t5*t6*t19*t22*wheelbase-V*fwheel_radius*phid*rwheel_radius*t2*t10*t12*t22*t25*wheelbase+com_to_rear*fwheel_radius*rwheel_radius*t4*t5*t9*t12*t13*t14*t27+com_to_rear*fwheel_radius*rwheel_radius*t2*t7*t8*t12*t20*t22*wheelbase+V*deltad*fwheel_radius*rwheel_radius*t2*t5*t9*t12*t14*t19*t22+V*deltad*fwheel_radius*rwheel_radius*t3*t4*t8*t11*t12*t26*wheelbase+V*deltad*fwheel_radius*rwheel_radius*t3*t4*t10*t12*t24*t26*wheelbase+V*fwheel_radius*phid*rwheel_radius*t3*t4*t9*t12*t14*t16*t24-V*fwheel_radius*phid*rwheel_radius*t2*t5*t10*t12*t19*t22*wheelbase*2.0+com_to_rear*fwheel_radius*g*rwheel_radius*t4*t5*t8*t12*t14*t20*t26+Ixz*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t4*t20*wheelbase+Ixz*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t20*t23*wheelbase-Ixx*Iyy*V*fwheel_radius*phid*rwheel_radius*t2*t5*t19*t22*wheelbase+Ixx*Izz*V*fwheel_radius*phid*rwheel_radius*t2*t5*t19*t22*wheelbase+Ixx*V*com_height*deltad*fwheel_radius*m*rwheel_radius*t5*t14*t17*t19+Ixx*V*com_height*fwheel_radius*m*phid*rwheel_radius*t4*t14*t16*t20*2.0+Ixx*V*deltad*fwheel_radius*m*rwheel_radius*t3*t4*t11*t26*wheelbase-Ixx*V*deltad*fwheel_radius*m*rwheel_radius*t3*t8*t23*t24*wheelbase+Iyy*V*deltad*fwheel_radius*m*rwheel_radius*t3*t4*t8*t26*wheelbase+Iyy*V*deltad*fwheel_radius*m*rwheel_radius*t4*t8*t16*t20*wheelbase+Izz*V*deltad*fwheel_radius*m*rwheel_radius*t4*t8*t20*t26*wheelbase+Ixx*com_height*com_to_rear*fwheel_radius*m*rwheel_radius*t4*t5*t13*t14*t27-Ixx*V*fwheel_radius*m*phid*rwheel_radius*t2*t8*t22*t25*wheelbase-Iyy*V*fwheel_radius*m*phid*rwheel_radius*t2*t5*t8*t22*wheelbase+Ixz*com_height*fwheel_radius*g*m*rwheel_radius*t4*t5*t14*t20*t26+Ixz*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t4*t16*t20*wheelbase+Ixx*V*com_height*deltad*fwheel_radius*m*rwheel_radius*t2*t5*t14*t19*t22+Ixx*V*com_height*fwheel_radius*m*phid*rwheel_radius*t3*t4*t14*t16*t24+Ixx*V*deltad*fwheel_radius*m*rwheel_radius*t3*t4*t8*t24*t26*wheelbase-Ixx*V*fwheel_radius*m*phid*rwheel_radius*t2*t5*t8*t19*t22*wheelbase*2.0-Iyy*V*fwheel_radius*m*phid*rwheel_radius*t2*t5*t8*t19*t22*wheelbase+Izz*V*fwheel_radius*m*phid*rwheel_radius*t2*t5*t8*t19*t22*wheelbase+V*fwheel_radius*phid*rwheel_radius*t2*t5*t8*t11*t12*t19*t22*wheelbase+Ixz*V*com_height*com_to_rear*fwheel_radius*m*phid*rwheel_radius*t2*t5*t19*t22*wheelbase*2.0)/(t2*t3*wheelbase*(Ixx*Iyy*fwheel_radius*rwheel_radius*t22+Ixx*fwheel_radius*m*rwheel_radius*t11+fwheel_radius*rwheel_radius*t8*t11*t12-fwheel_radius*rwheel_radius*t6*t19*t22+fwheel_radius*rwheel_radius*t8*t11*t12*t26+fwheel_radius*rwheel_radius*t8*t12*t14*t27+fwheel_radius*rwheel_radius*t10*t12*t22*t24-Ixx*Iyy*fwheel_radius*rwheel_radius*t19*t22+Ixx*Izz*fwheel_radius*rwheel_radius*t19*t22+Ixx*fwheel_radius*m*rwheel_radius*t11*t26+Ixx*fwheel_radius*m*rwheel_radius*t14*t27+Iyy*fwheel_radius*m*rwheel_radius*t8*t22+Ixx*fwheel_radius*m*rwheel_radius*t8*t22*t24-Iyy*fwheel_radius*m*rwheel_radius*t8*t19*t22+Izz*fwheel_radius*m*rwheel_radius*t8*t19*t22-fwheel_radius*rwheel_radius*t8*t11*t12*t19*t22-Ixz*com_height*com_to_rear*fwheel_radius*m*rwheel_radius*t19*t22*2.0-fwheel_radius*rwheel_radius*t2*t3*t4*t5*t9*t12*wheelbase*2.0-Ixx*com_height*fwheel_radius*m*rwheel_radius*t2*t3*t4*t5*wheelbase*2.0));
}

adouble phidd_eqn(adouble& phi, adouble& V, adouble& phid, adouble& delta, adouble& deltad, adouble& torque_front, adouble& torque_rear) {
  /*
   state -- state vector, order [phi, V, phid, delta]
% controls -- control vector, order [torque_front torque_rear, deltad]
% sys_params -- parameter vector, order [g, m, Ixx, Iyy, Izz, Ixz, wheelbase, com_to_rear, com_height, fwheel_radius, rwheel_radius]
  */
  adouble t2 = cos(delta);
  adouble t3 = cos(phi);
  adouble t4 = sin(delta);
  adouble t5 = sin(phi);
  adouble t6 = Ixz*Ixz;
  adouble t7 = Iyy*Iyy;
  adouble t8 = V*V;
  adouble t9 = com_height*com_height;
  adouble t10 = com_height*com_height*com_height;
  adouble t12 = com_to_rear*com_to_rear;
  adouble t13 = com_to_rear*com_to_rear*com_to_rear;
  adouble t14 = m*m;
  adouble t15 = phid*phid;
  adouble t16 = wheelbase*wheelbase;
  adouble t17 = wheelbase*wheelbase*wheelbase;
  adouble t11 = t9*t9;
  adouble t18 = t2*t2;
  adouble t19 = t2*t2*t2;
  adouble t21 = t3*t3;
  adouble t22 = t3*t3*t3;
  adouble t24 = t4*t4;
  adouble t25 = t4*t4*t4;
  adouble t27 = t5*t5;
  adouble t28 = t5*t5*t5;
  adouble t20 = t18*t18;
  adouble t23 = t21*t21;
  adouble t26 = t24*t24;
  adouble t29 = -t18;
  adouble t30 = t18*t21;
  adouble t31 = t29+t30+1.0;
  adouble t32 = 1.0/sqrt(t31);
  return -(-fwheel_radius*rwheel_radius*t5*t7*t8*t26+Iyy*Izz*fwheel_radius*rwheel_radius*t5*t8*t26+Ixz*fwheel_radius*t4*t17*t19*t22*torque_rear+Ixz*rwheel_radius*t4*t17*t30*t32*torque_front+fwheel_radius*rwheel_radius*t5*t7*t8*t21*t26-fwheel_radius*rwheel_radius*t8*t11*t14*t26*t28+(Izz*Izz)*fwheel_radius*rwheel_radius*t5*t8*t21*t26-Iyy*Izz*fwheel_radius*rwheel_radius*t5*t8*t21*t26*2.0-Iyy*fwheel_radius*m*rwheel_radius*t5*t8*t9*t26-Iyy*fwheel_radius*m*rwheel_radius*t5*t8*t12*t24-Iyy*fwheel_radius*m*rwheel_radius*t8*t9*t26*t28+Izz*fwheel_radius*m*rwheel_radius*t5*t8*t12*t24+Izz*fwheel_radius*m*rwheel_radius*t8*t9*t26*t28-Ixz*rwheel_radius*t4*t17*t20*t21*t32*torque_front+Ixz*rwheel_radius*t4*t17*t20*t23*t32*torque_front-fwheel_radius*rwheel_radius*t5*t8*t9*t12*t14*t24+fwheel_radius*rwheel_radius*t5*t6*t15*t16*t24*t30+Ixz*V*deltad*fwheel_radius*m*rwheel_radius*t17*t20*t22+V*com_height*deltad*fwheel_radius*rwheel_radius*t3*t13*t14*wheelbase+Iyy*fwheel_radius*m*rwheel_radius*t5*t8*t12*t18*t24+Iyy*fwheel_radius*m*rwheel_radius*t5*t8*t9*t21*t26-Izz*fwheel_radius*m*rwheel_radius*t5*t8*t9*t21*t26+Izz*fwheel_radius*m*rwheel_radius*t5*t8*t12*t24*t29+Izz*fwheel_radius*m*rwheel_radius*t5*t8*t16*t24*t30+com_height*com_to_rear*fwheel_radius*m*t4*t17*t19*t22*torque_rear+com_height*com_to_rear*m*rwheel_radius*t4*t17*t30*t32*torque_front+com_height*fwheel_radius*rwheel_radius*t4*t8*t14*t17*t19*t22+fwheel_radius*rwheel_radius*t5*t8*t9*t12*t14*t18*t24-fwheel_radius*rwheel_radius*t5*t8*t9*t14*t16*t24*t30*3.0-com_height*fwheel_radius*g*rwheel_radius*t5*t14*(t16*t16)*t20*t22+Ixz*V*deltad*fwheel_radius*m*rwheel_radius*t17*t18*t22*t24+Ixz*V*deltad*fwheel_radius*m*rwheel_radius*t3*t12*t18*wheelbase-Ixz*V*deltad*fwheel_radius*m*rwheel_radius*t3*t12*t20*wheelbase+Ixz*V*deltad*fwheel_radius*m*rwheel_radius*t3*t12*t24*wheelbase-Ixz*V*deltad*fwheel_radius*m*rwheel_radius*t3*t12*t26*wheelbase+V*com_height*com_to_rear*deltad*fwheel_radius*rwheel_radius*t14*t17*t18*t22+Iyy*com_height*fwheel_radius*m*rwheel_radius*t2*t3*t8*t25*wheelbase-Iyy*com_height*fwheel_radius*m*rwheel_radius*t2*t8*t22*t25*wheelbase+Izz*com_height*fwheel_radius*m*rwheel_radius*t2*t8*t22*t25*wheelbase-V*com_height*deltad*fwheel_radius*rwheel_radius*t3*t13*t14*t26*wheelbase+V*com_height*deltad*fwheel_radius*rwheel_radius*t3*t13*t14*t29*wheelbase+Iyy*fwheel_radius*m*rwheel_radius*t5*t8*t16*t21*t24*t29+com_height*fwheel_radius*g*rwheel_radius*t3*t5*t12*t14*t16*t20+com_height*fwheel_radius*g*rwheel_radius*t3*t5*t12*t14*t16*t29-com_height*com_to_rear*m*rwheel_radius*t4*t17*t20*t21*t32*torque_front+com_height*com_to_rear*m*rwheel_radius*t4*t17*t20*t23*t32*torque_front+com_height*fwheel_radius*rwheel_radius*t2*t3*t4*t8*t12*t14*wheelbase-com_height*fwheel_radius*rwheel_radius*t3*t4*t8*t12*t14*t19*wheelbase+fwheel_radius*g*rwheel_radius*t4*t9*t14*t17*t19*t21*t27*2.0+fwheel_radius*g*rwheel_radius*t3*t10*t14*t16*t24*t28*t29+fwheel_radius*rwheel_radius*t5*t9*t12*t14*t15*t16*t24*t30+fwheel_radius*rwheel_radius*t2*t3*t8*t10*t14*t25*t27*wheelbase*3.0+Iyy*com_height*fwheel_radius*m*rwheel_radius*t2*t3*t8*t25*t27*wheelbase*2.0-Izz*com_height*fwheel_radius*m*rwheel_radius*t2*t3*t8*t25*t27*wheelbase*2.0+V*com_height*deltad*fwheel_radius*rwheel_radius*t3*t13*t14*t24*t29*wheelbase+V*com_to_rear*deltad*fwheel_radius*rwheel_radius*t3*t10*t14*t24*t27*wheelbase-V*com_to_rear*deltad*fwheel_radius*rwheel_radius*t3*t10*t14*t26*t27*wheelbase+V*com_to_rear*fwheel_radius*phid*rwheel_radius*t9*t14*t16*t18*t22*t24*2.0+Iyy*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t3*t24*wheelbase-Iyy*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t3*t26*wheelbase-Iyy*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t22*t24*wheelbase+Iyy*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t22*t26*wheelbase+Izz*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t22*t24*wheelbase-Izz*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t22*t26*wheelbase-Ixz*Iyy*V*fwheel_radius*phid*rwheel_radius*t2*t5*t21*t25*wheelbase*2.0+Ixz*Izz*V*fwheel_radius*phid*rwheel_radius*t2*t5*t21*t25*wheelbase*2.0+Ixz*V*com_height*fwheel_radius*m*phid*rwheel_radius*t16*t18*t22*t24*2.0-Ixz*V*deltad*fwheel_radius*m*rwheel_radius*t3*t12*t18*t24*wheelbase*2.0+Ixz*V*fwheel_radius*m*phid*rwheel_radius*t4*t5*t17*t19*t21+Ixz*com_height*com_to_rear*fwheel_radius*m*rwheel_radius*t5*t15*t16*t24*t30*2.0+Iyy*com_height*fwheel_radius*g*m*rwheel_radius*t3*t5*t16*t24*t29+Iyy*com_height*fwheel_radius*g*m*rwheel_radius*t5*t16*t18*t22*t24+Izz*com_height*fwheel_radius*g*m*rwheel_radius*t5*t16*t22*t24*t29+Iyy*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t3*t24*t29*wheelbase+Iyy*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t18*t22*t24*wheelbase+Izz*V*com_height*com_to_rear*deltad*fwheel_radius*m*rwheel_radius*t22*t24*t29*wheelbase-Ixz*V*com_height*deltad*fwheel_radius*m*rwheel_radius*t4*t5*t16*t19*t21-Ixz*V*com_height*deltad*fwheel_radius*m*rwheel_radius*t2*t5*t16*t21*t25+Ixz*V*com_height*fwheel_radius*m*phid*rwheel_radius*t3*t16*t24*t27*t29-Ixz*V*fwheel_radius*m*phid*rwheel_radius*t2*t5*t9*t21*t25*wheelbase*2.0+V*com_height*com_to_rear*fwheel_radius*phid*rwheel_radius*t4*t5*t14*t17*t19*t21-V*com_to_rear*deltad*fwheel_radius*rwheel_radius*t2*t4*t5*t9*t14*t16*t21*2.0+V*com_to_rear*deltad*fwheel_radius*rwheel_radius*t4*t5*t9*t14*t16*t19*t21+V*com_to_rear*deltad*fwheel_radius*rwheel_radius*t2*t5*t9*t14*t16*t21*t25+V*com_to_rear*deltad*fwheel_radius*rwheel_radius*t3*t10*t14*t24*t27*t29*wheelbase+V*com_to_rear*fwheel_radius*phid*rwheel_radius*t3*t9*t14*t16*t24*t27*t29-V*com_to_rear*fwheel_radius*phid*rwheel_radius*t2*t5*t10*t14*t21*t25*wheelbase*2.0-Iyy*V*com_height*com_to_rear*fwheel_radius*m*phid*rwheel_radius*t2*t5*t21*t25*wheelbase*2.0+Izz*V*com_height*com_to_rear*fwheel_radius*m*phid*rwheel_radius*t2*t5*t21*t25*wheelbase*2.0)/(wheelbase*(Ixx*fwheel_radius*m*rwheel_radius*t17*t20*t22+fwheel_radius*rwheel_radius*t9*t14*t17*t20*t22+fwheel_radius*rwheel_radius*t6*t22*t24*t29*wheelbase+Ixx*Iyy*fwheel_radius*rwheel_radius*t3*t18*t24*wheelbase+Ixx*Iyy*fwheel_radius*rwheel_radius*t22*t24*t29*wheelbase+Ixx*Izz*fwheel_radius*rwheel_radius*t18*t22*t24*wheelbase+Ixx*fwheel_radius*m*rwheel_radius*t3*t12*t18*wheelbase-Ixx*fwheel_radius*m*rwheel_radius*t3*t12*t20*wheelbase+fwheel_radius*rwheel_radius*t3*t9*t12*t14*t18*wheelbase-fwheel_radius*rwheel_radius*t3*t9*t12*t14*t20*wheelbase+Iyy*fwheel_radius*m*rwheel_radius*t3*t9*t18*t24*wheelbase+Iyy*fwheel_radius*m*rwheel_radius*t9*t22*t24*t29*wheelbase+Izz*fwheel_radius*m*rwheel_radius*t9*t18*t22*t24*wheelbase-fwheel_radius*rwheel_radius*t4*t5*t10*t14*t16*t19*t21*2.0+fwheel_radius*rwheel_radius*t3*t11*t14*t18*t24*t27*wheelbase+fwheel_radius*rwheel_radius*t9*t12*t14*t22*t24*t29*wheelbase-Ixz*com_height*com_to_rear*fwheel_radius*m*rwheel_radius*t18*t22*t24*wheelbase*2.0-Ixx*com_height*fwheel_radius*m*rwheel_radius*t4*t5*t16*t19*t21*2.0+Ixx*fwheel_radius*m*rwheel_radius*t3*t9*t18*t24*t27*wheelbase));

}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    return states[4]*0. ;//controls[0] * controls[0] + controls[1] * controls[1] + controls[2] * controls[2]; // Squared torque front, rear, and steer acceleration.
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble xdot, ydot, vdot;

   adouble phi = states[0];
   adouble V = states[1];
   adouble phid = states[2];
   adouble delta = states[3];
   adouble deltad = states[4];
   adouble psi = states[5];
   adouble xr = states[6];
   adouble yr = states[7];

   adouble tr = controls[0];
   adouble tf = controls[1];
   adouble deltadd = controls[2];
   
   derivatives[0] = phid;
   derivatives[1] = Vd_eqn(phi, V, phid, delta,  deltad,  tf,  tr); 
   // cout << phi << ", " << V << ", " << phid << ", " << delta << ", " << deltad << ", " << tf << ", " << tr << "\n";
   derivatives[2] = phidd_eqn(phi, V, phid, delta,  deltad,  tf,  tr); 
   derivatives[3] = deltad;
   derivatives[4] = deltadd;
   derivatives[5] = psid_eqn(phi, V, delta);
   derivatives[6] = V * cos(psi);
   // cout << derivatives[5] << "\n";
   derivatives[7] = V * sin(psi);
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble x10 = initial_states[ 0 ];
   adouble x20 = initial_states[ 1 ];
   adouble x30 = initial_states[ 2 ];
   adouble x40 = initial_states[ 3 ];
  
   adouble x50 = initial_states[ 4 ];
   adouble x60 = initial_states[ 5 ];
   adouble x70 = initial_states[ 6 ];
   adouble x80 = initial_states[ 7 ];

   adouble x1f = final_states[ 0 ];
   adouble x2f = final_states[ 1 ];
   adouble x3f = final_states[ 2 ];
   adouble x4f = final_states[ 3 ];

   adouble x5f = final_states[ 4 ];
   adouble x6f = final_states[ 5 ];
   adouble x7f = final_states[ 6 ];
   adouble x8f = final_states[ 7 ];

   e[ 0 ] = x10;
   e[ 1 ] = x20;
   e[ 2 ] = x30;
   e[ 3 ] = x40;
   e[ 4 ] = x50;
   e[ 5 ] = x60;
   e[ 6 ] = x70;
   e[ 7 ] = x80;
  
   e[ 8 ] = x1f;
   e[ 9 ] = x2f;
   e[ 10 ] = x3f;
   e[ 11 ] = x4f;
   e[ 12 ] = x5f;
   e[ 13 ] = x6f;
   e[ 14 ] = x7f;
   e[ 15 ] = x8f;
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////


int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        						= "bike";

    problem.outfilename                 	= "bike.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   							= 1;
    problem.nlinkages                   	= 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   				= 8;
    problem.phases(1).ncontrols 				= 3;
    problem.phases(1).nevents   				= 16;
    problem.phases(1).npath     				= 0;
    int nodes = 10;
    problem.phases(1).nodes           << nodes;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.states(0) = -1.57; // phi
    problem.phases(1).bounds.lower.states(1) = 0.0; // V
    problem.phases(1).bounds.lower.states(2) = -30.0; // phid
    problem.phases(1).bounds.lower.states(3) = -1.57; // delta
    problem.phases(1).bounds.lower.states(4) = -30.0; // deltad
    problem.phases(1).bounds.lower.states(5) = -3.142; // psi
    problem.phases(1).bounds.lower.states(6) = -20.0; // xr
    problem.phases(1).bounds.lower.states(7) = -20.0; // yr

    problem.phases(1).bounds.upper.states(0) = 1.57;
    problem.phases(1).bounds.upper.states(1) = 2.0;
    problem.phases(1).bounds.upper.states(2) = 30.0;
    problem.phases(1).bounds.upper.states(3) = 1.57;
    problem.phases(1).bounds.upper.states(4) = 30.0;
    problem.phases(1).bounds.upper.states(5) = 3.142;
    problem.phases(1).bounds.upper.states(6) = 20.0;
    problem.phases(1).bounds.upper.states(7) = 20.0;

    problem.phases(1).bounds.lower.controls(0) = -20.0; // tfront
    problem.phases(1).bounds.lower.controls(1) = -20.0; // trear
    problem.phases(1).bounds.lower.controls(2) = -5.; // deltadd

    problem.phases(1).bounds.upper.controls(0) = 20.0;
    problem.phases(1).bounds.upper.controls(1) = 20.0;
    problem.phases(1).bounds.upper.controls(2) = 5.0;

    // Initial.
    problem.phases(1).bounds.lower.events(0) = 0.05;
    problem.phases(1).bounds.lower.events(1) = 2.0;
    problem.phases(1).bounds.lower.events(2) = 0.0;
    problem.phases(1).bounds.lower.events(3) = 0.0;
    problem.phases(1).bounds.lower.events(4) = 0.0;
    problem.phases(1).bounds.lower.events(5) = 0.0;
    problem.phases(1).bounds.lower.events(6) = 0.0;
    problem.phases(1).bounds.lower.events(7) = 0.0;

    // Final.
    problem.phases(1).bounds.lower.events(8) = 0.0;
    problem.phases(1).bounds.lower.events(9) = 0.0;
    problem.phases(1).bounds.lower.events(10) = 0.0;
    problem.phases(1).bounds.lower.events(11) = 1.0;
    problem.phases(1).bounds.lower.events(12) = 0.0;
    problem.phases(1).bounds.lower.events(13) = 0.0;
    problem.phases(1).bounds.lower.events(14) = 4.0;
    problem.phases(1).bounds.lower.events(15) = 0.0;

    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 3.;
    problem.phases(1).bounds.upper.EndTime      = 6.;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 				= &integrand_cost;
    problem.endpoint_cost 					= &endpoint_cost;
    problem.dae 								= &dae;
    problem.events 							= &events;
    problem.linkages							= &linkages;



////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd x0(8, nodes);

    x0 <<  linspace(0.0, 0.0, nodes), // phi
           linspace(2.0, 0.0, nodes), // V
           linspace(0.0, 0.0, nodes), // phid
           linspace(0.0, 1.0, nodes), // delta
           linspace(0.0, 0.0, nodes), // deltad
           linspace(0.0, 0.0, nodes), // yaw
           linspace(0.0, 10.0, nodes), // xr
           linspace(0.0, 0.0, nodes); // yr

    problem.phases(1).guess.controls       = zeros(3, nodes);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 3.0, nodes);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.collocation_method = "Hermite-Simpson";
////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x 		= solution.get_states_in_phase(1);
    u 		= solution.get_controls_in_phase(1);
    t 		= solution.get_time_in_phase(1);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4");

    plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2");


    plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4",
                                  "pdf", "twolinkarm_states.pdf");

    plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2",
                              "pdf", "twolinkarm_controls.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

