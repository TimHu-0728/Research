function [duxva,duyva,duxxva,duyyva] = collocation(in1,in2,in3,in4,in5)
%COLLOCATION
%    [DUXVA,DUYVA,DUXXVA,DUYYVA] = COLLOCATION(IN1,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    20-Mar-2017 12:05:48

x1r1 = in1(:,1);
x1r2 = in1(:,2);
x1r3 = in1(:,3);
x1r4 = in1(:,4);
x1r5 = in1(:,5);
x2r1 = in3(:,1);
x2r2 = in3(:,2);
x2r3 = in3(:,3);
x2r4 = in3(:,4);
x2r5 = in3(:,5);
x1ry1r1 = in5(:,1);
x1ry1r2 = in5(:,2);
x1ry1r3 = in5(:,3);
x1ry1r4 = in5(:,4);
x1ry1r5 = in5(:,5);
y1r1 = in2(:,1);
y1r2 = in2(:,2);
y1r3 = in2(:,3);
y1r4 = in2(:,4);
y1r5 = in2(:,5);
y2r1 = in4(:,1);
y2r2 = in4(:,2);
y2r3 = in4(:,3);
y2r4 = in4(:,4);
y2r5 = in4(:,5);
t2 = x2r2.*x1ry1r3.*y1r4.*y2r5;
t3 = x2r2.*x1ry1r4.*y1r5.*y2r3;
t4 = x2r2.*x1ry1r5.*y1r3.*y2r4;
t5 = x2r3.*x1ry1r2.*y1r5.*y2r4;
t6 = x2r3.*x1ry1r4.*y1r2.*y2r5;
t7 = x2r3.*x1ry1r5.*y1r4.*y2r2;
t8 = x2r4.*x1ry1r2.*y1r3.*y2r5;
t9 = x2r4.*x1ry1r3.*y1r5.*y2r2;
t10 = x2r4.*x1ry1r5.*y1r2.*y2r3;
t11 = x2r5.*x1ry1r2.*y1r4.*y2r3;
t12 = x2r5.*x1ry1r3.*y1r2.*y2r4;
t13 = x2r5.*x1ry1r4.*y1r3.*y2r2;
t14 = x1r1.*x2r2.*x1ry1r3.*y1r4.*y2r5;
t15 = x1r1.*x2r2.*x1ry1r4.*y1r5.*y2r3;
t16 = x1r1.*x2r2.*x1ry1r5.*y1r3.*y2r4;
t17 = x1r1.*x2r3.*x1ry1r2.*y1r5.*y2r4;
t18 = x1r1.*x2r3.*x1ry1r4.*y1r2.*y2r5;
t19 = x1r1.*x2r3.*x1ry1r5.*y1r4.*y2r2;
t20 = x1r1.*x2r4.*x1ry1r2.*y1r3.*y2r5;
t21 = x1r1.*x2r4.*x1ry1r3.*y1r5.*y2r2;
t22 = x1r1.*x2r4.*x1ry1r5.*y1r2.*y2r3;
t23 = x1r1.*x2r5.*x1ry1r2.*y1r4.*y2r3;
t24 = x1r1.*x2r5.*x1ry1r3.*y1r2.*y2r4;
t25 = x1r1.*x2r5.*x1ry1r4.*y1r3.*y2r2;
t26 = x1r2.*x2r1.*x1ry1r3.*y1r5.*y2r4;
t27 = x1r2.*x2r1.*x1ry1r4.*y1r3.*y2r5;
t28 = x1r2.*x2r1.*x1ry1r5.*y1r4.*y2r3;
t29 = x1r2.*x2r3.*x1ry1r1.*y1r4.*y2r5;
t30 = x1r2.*x2r3.*x1ry1r4.*y1r5.*y2r1;
t31 = x1r2.*x2r3.*x1ry1r5.*y1r1.*y2r4;
t32 = x1r2.*x2r4.*x1ry1r1.*y1r5.*y2r3;
t33 = x1r2.*x2r4.*x1ry1r3.*y1r1.*y2r5;
t34 = x1r2.*x2r4.*x1ry1r5.*y1r3.*y2r1;
t35 = x1r2.*x2r5.*x1ry1r1.*y1r3.*y2r4;
t36 = x1r2.*x2r5.*x1ry1r3.*y1r4.*y2r1;
t37 = x1r2.*x2r5.*x1ry1r4.*y1r1.*y2r3;
t38 = x1r3.*x2r1.*x1ry1r2.*y1r4.*y2r5;
t39 = x1r3.*x2r1.*x1ry1r4.*y1r5.*y2r2;
t40 = x1r3.*x2r1.*x1ry1r5.*y1r2.*y2r4;
t41 = x1r3.*x2r2.*x1ry1r1.*y1r5.*y2r4;
t42 = x1r3.*x2r2.*x1ry1r4.*y1r1.*y2r5;
t43 = x1r3.*x2r2.*x1ry1r5.*y1r4.*y2r1;
t44 = x1r3.*x2r4.*x1ry1r1.*y1r2.*y2r5;
t45 = x1r3.*x2r4.*x1ry1r2.*y1r5.*y2r1;
t46 = x1r3.*x2r4.*x1ry1r5.*y1r1.*y2r2;
t47 = x1r3.*x2r5.*x1ry1r1.*y1r4.*y2r2;
t48 = x1r3.*x2r5.*x1ry1r2.*y1r1.*y2r4;
t49 = x1r3.*x2r5.*x1ry1r4.*y1r2.*y2r1;
t50 = x1r4.*x2r1.*x1ry1r2.*y1r5.*y2r3;
t51 = x1r4.*x2r1.*x1ry1r3.*y1r2.*y2r5;
t52 = x1r4.*x2r1.*x1ry1r5.*y1r3.*y2r2;
t53 = x1r4.*x2r2.*x1ry1r1.*y1r3.*y2r5;
t54 = x1r4.*x2r2.*x1ry1r3.*y1r5.*y2r1;
t55 = x1r4.*x2r2.*x1ry1r5.*y1r1.*y2r3;
t56 = x1r4.*x2r3.*x1ry1r1.*y1r5.*y2r2;
t57 = x1r4.*x2r3.*x1ry1r2.*y1r1.*y2r5;
t58 = x1r4.*x2r3.*x1ry1r5.*y1r2.*y2r1;
t59 = x1r4.*x2r5.*x1ry1r1.*y1r2.*y2r3;
t60 = x1r4.*x2r5.*x1ry1r2.*y1r3.*y2r1;
t61 = x1r4.*x2r5.*x1ry1r3.*y1r1.*y2r2;
t62 = x1r5.*x2r1.*x1ry1r2.*y1r3.*y2r4;
t63 = x1r5.*x2r1.*x1ry1r3.*y1r4.*y2r2;
t64 = x1r5.*x2r1.*x1ry1r4.*y1r2.*y2r3;
t65 = x1r5.*x2r2.*x1ry1r1.*y1r4.*y2r3;
t66 = x1r5.*x2r2.*x1ry1r3.*y1r1.*y2r4;
t67 = x1r5.*x2r2.*x1ry1r4.*y1r3.*y2r1;
t68 = x1r5.*x2r3.*x1ry1r1.*y1r2.*y2r4;
t69 = x1r5.*x2r3.*x1ry1r2.*y1r4.*y2r1;
t70 = x1r5.*x2r3.*x1ry1r4.*y1r1.*y2r2;
t71 = x1r5.*x2r4.*x1ry1r1.*y1r3.*y2r2;
t72 = x1r5.*x2r4.*x1ry1r2.*y1r1.*y2r3;
t73 = x1r5.*x2r4.*x1ry1r3.*y1r2.*y2r1;
t88 = x1r1.*x2r2.*x1ry1r3.*y1r5.*y2r4;
t89 = x1r1.*x2r2.*x1ry1r4.*y1r3.*y2r5;
t90 = x1r1.*x2r2.*x1ry1r5.*y1r4.*y2r3;
t91 = x1r1.*x2r3.*x1ry1r2.*y1r4.*y2r5;
t92 = x1r1.*x2r3.*x1ry1r4.*y1r5.*y2r2;
t93 = x1r1.*x2r3.*x1ry1r5.*y1r2.*y2r4;
t94 = x1r1.*x2r4.*x1ry1r2.*y1r5.*y2r3;
t95 = x1r1.*x2r4.*x1ry1r3.*y1r2.*y2r5;
t96 = x1r1.*x2r4.*x1ry1r5.*y1r3.*y2r2;
t97 = x1r1.*x2r5.*x1ry1r2.*y1r3.*y2r4;
t98 = x1r1.*x2r5.*x1ry1r3.*y1r4.*y2r2;
t99 = x1r1.*x2r5.*x1ry1r4.*y1r2.*y2r3;
t100 = x1r2.*x2r1.*x1ry1r3.*y1r4.*y2r5;
t101 = x1r2.*x2r1.*x1ry1r4.*y1r5.*y2r3;
t102 = x1r2.*x2r1.*x1ry1r5.*y1r3.*y2r4;
t103 = x1r2.*x2r3.*x1ry1r1.*y1r5.*y2r4;
t104 = x1r2.*x2r3.*x1ry1r4.*y1r1.*y2r5;
t105 = x1r2.*x2r3.*x1ry1r5.*y1r4.*y2r1;
t106 = x1r2.*x2r4.*x1ry1r1.*y1r3.*y2r5;
t107 = x1r2.*x2r4.*x1ry1r3.*y1r5.*y2r1;
t108 = x1r2.*x2r4.*x1ry1r5.*y1r1.*y2r3;
t109 = x1r2.*x2r5.*x1ry1r1.*y1r4.*y2r3;
t110 = x1r2.*x2r5.*x1ry1r3.*y1r1.*y2r4;
t111 = x1r2.*x2r5.*x1ry1r4.*y1r3.*y2r1;
t112 = x1r3.*x2r1.*x1ry1r2.*y1r5.*y2r4;
t113 = x1r3.*x2r1.*x1ry1r4.*y1r2.*y2r5;
t114 = x1r3.*x2r1.*x1ry1r5.*y1r4.*y2r2;
t115 = x1r3.*x2r2.*x1ry1r1.*y1r4.*y2r5;
t116 = x1r3.*x2r2.*x1ry1r4.*y1r5.*y2r1;
t117 = x1r3.*x2r2.*x1ry1r5.*y1r1.*y2r4;
t118 = x1r3.*x2r4.*x1ry1r1.*y1r5.*y2r2;
t119 = x1r3.*x2r4.*x1ry1r2.*y1r1.*y2r5;
t120 = x1r3.*x2r4.*x1ry1r5.*y1r2.*y2r1;
t121 = x1r3.*x2r5.*x1ry1r1.*y1r2.*y2r4;
t122 = x1r3.*x2r5.*x1ry1r2.*y1r4.*y2r1;
t123 = x1r3.*x2r5.*x1ry1r4.*y1r1.*y2r2;
t124 = x1r4.*x2r1.*x1ry1r2.*y1r3.*y2r5;
t125 = x1r4.*x2r1.*x1ry1r3.*y1r5.*y2r2;
t126 = x1r4.*x2r1.*x1ry1r5.*y1r2.*y2r3;
t127 = x1r4.*x2r2.*x1ry1r1.*y1r5.*y2r3;
t128 = x1r4.*x2r2.*x1ry1r3.*y1r1.*y2r5;
t129 = x1r4.*x2r2.*x1ry1r5.*y1r3.*y2r1;
t130 = x1r4.*x2r3.*x1ry1r1.*y1r2.*y2r5;
t131 = x1r4.*x2r3.*x1ry1r2.*y1r5.*y2r1;
t132 = x1r4.*x2r3.*x1ry1r5.*y1r1.*y2r2;
t133 = x1r4.*x2r5.*x1ry1r1.*y1r3.*y2r2;
t134 = x1r4.*x2r5.*x1ry1r2.*y1r1.*y2r3;
t135 = x1r4.*x2r5.*x1ry1r3.*y1r2.*y2r1;
t136 = x1r5.*x2r1.*x1ry1r2.*y1r4.*y2r3;
t137 = x1r5.*x2r1.*x1ry1r3.*y1r2.*y2r4;
t138 = x1r5.*x2r1.*x1ry1r4.*y1r3.*y2r2;
t139 = x1r5.*x2r2.*x1ry1r1.*y1r3.*y2r4;
t140 = x1r5.*x2r2.*x1ry1r3.*y1r4.*y2r1;
t141 = x1r5.*x2r2.*x1ry1r4.*y1r1.*y2r3;
t142 = x1r5.*x2r3.*x1ry1r1.*y1r4.*y2r2;
t143 = x1r5.*x2r3.*x1ry1r2.*y1r1.*y2r4;
t144 = x1r5.*x2r3.*x1ry1r4.*y1r2.*y2r1;
t145 = x1r5.*x2r4.*x1ry1r1.*y1r2.*y2r3;
t146 = x1r5.*x2r4.*x1ry1r2.*y1r3.*y2r1;
t147 = x1r5.*x2r4.*x1ry1r3.*y1r1.*y2r2;
t74 = t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73-t88-t89-t90-t91-t92-t93-t94-t95-t96-t97-t98-t99-t100-t101-t102-t103-t104-t105-t106-t107-t108-t109-t110-t111-t112-t113-t114-t115-t116-t117-t118-t119-t120-t121-t122-t123-t124-t125-t126-t127-t128-t129-t130-t131-t132-t133-t134-t135-t136-t137-t138-t139-t140-t141-t142-t143-t144-t145-t146-t147;
t75 = 1.0./t74;
t76 = x2r1.*x1ry1r3.*y1r5.*y2r4;
t77 = x2r1.*x1ry1r4.*y1r3.*y2r5;
t78 = x2r1.*x1ry1r5.*y1r4.*y2r3;
t79 = x2r3.*x1ry1r1.*y1r4.*y2r5;
t80 = x2r3.*x1ry1r4.*y1r5.*y2r1;
t81 = x2r3.*x1ry1r5.*y1r1.*y2r4;
t82 = x2r4.*x1ry1r1.*y1r5.*y2r3;
t83 = x2r4.*x1ry1r3.*y1r1.*y2r5;
t84 = x2r4.*x1ry1r5.*y1r3.*y2r1;
t85 = x2r5.*x1ry1r1.*y1r3.*y2r4;
t86 = x2r5.*x1ry1r3.*y1r4.*y2r1;
t87 = x2r5.*x1ry1r4.*y1r1.*y2r3;
t148 = x2r1.*x1ry1r2.*y1r4.*y2r5;
t149 = x2r1.*x1ry1r4.*y1r5.*y2r2;
t150 = x2r1.*x1ry1r5.*y1r2.*y2r4;
t151 = x2r2.*x1ry1r1.*y1r5.*y2r4;
t152 = x2r2.*x1ry1r4.*y1r1.*y2r5;
t153 = x2r2.*x1ry1r5.*y1r4.*y2r1;
t154 = x2r4.*x1ry1r1.*y1r2.*y2r5;
t155 = x2r4.*x1ry1r2.*y1r5.*y2r1;
t156 = x2r4.*x1ry1r5.*y1r1.*y2r2;
t157 = x2r5.*x1ry1r1.*y1r4.*y2r2;
t158 = x2r5.*x1ry1r2.*y1r1.*y2r4;
t159 = x2r5.*x1ry1r4.*y1r2.*y2r1;
t160 = x2r1.*x1ry1r2.*y1r5.*y2r3;
t161 = x2r1.*x1ry1r3.*y1r2.*y2r5;
t162 = x2r1.*x1ry1r5.*y1r3.*y2r2;
t163 = x2r2.*x1ry1r1.*y1r3.*y2r5;
t164 = x2r2.*x1ry1r3.*y1r5.*y2r1;
t165 = x2r2.*x1ry1r5.*y1r1.*y2r3;
t166 = x2r3.*x1ry1r1.*y1r5.*y2r2;
t167 = x2r3.*x1ry1r2.*y1r1.*y2r5;
t168 = x2r3.*x1ry1r5.*y1r2.*y2r1;
t169 = x2r5.*x1ry1r1.*y1r2.*y2r3;
t170 = x2r5.*x1ry1r2.*y1r3.*y2r1;
t171 = x2r5.*x1ry1r3.*y1r1.*y2r2;
t172 = x2r1.*x1ry1r2.*y1r3.*y2r4;
t173 = x2r1.*x1ry1r3.*y1r4.*y2r2;
t174 = x2r1.*x1ry1r4.*y1r2.*y2r3;
t175 = x2r2.*x1ry1r1.*y1r4.*y2r3;
t176 = x2r2.*x1ry1r3.*y1r1.*y2r4;
t177 = x2r2.*x1ry1r4.*y1r3.*y2r1;
t178 = x2r3.*x1ry1r1.*y1r2.*y2r4;
t179 = x2r3.*x1ry1r2.*y1r4.*y2r1;
t180 = x2r3.*x1ry1r4.*y1r1.*y2r2;
t181 = x2r4.*x1ry1r1.*y1r3.*y2r2;
t182 = x2r4.*x1ry1r2.*y1r1.*y2r3;
t183 = x2r4.*x1ry1r3.*y1r2.*y2r1;
duxva = [-t75.*(t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t86+t87+t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162+t163+t164+t165+t166+t167+t168+t169+t170+t171+t172+t173+t174+t175+t176+t177+t178+t179+t180+t181+t182+t183-x2r1.*x1ry1r2.*y1r4.*y2r3-x2r1.*x1ry1r3.*y1r2.*y2r4-x2r1.*x1ry1r4.*y1r3.*y2r2-x2r2.*x1ry1r1.*y1r3.*y2r4-x2r2.*x1ry1r3.*y1r4.*y2r1-x2r2.*x1ry1r4.*y1r1.*y2r3-x2r3.*x1ry1r1.*y1r4.*y2r2-x2r3.*x1ry1r2.*y1r1.*y2r4-x2r3.*x1ry1r4.*y1r2.*y2r1-x2r4.*x1ry1r1.*y1r2.*y2r3-x2r4.*x1ry1r2.*y1r3.*y2r1-x2r4.*x1ry1r3.*y1r1.*y2r2-x2r1.*x1ry1r2.*y1r3.*y2r5-x2r1.*x1ry1r3.*y1r5.*y2r2-x2r1.*x1ry1r5.*y1r2.*y2r3-x2r2.*x1ry1r1.*y1r5.*y2r3-x2r2.*x1ry1r3.*y1r1.*y2r5-x2r2.*x1ry1r5.*y1r3.*y2r1-x2r3.*x1ry1r1.*y1r2.*y2r5-x2r3.*x1ry1r2.*y1r5.*y2r1-x2r3.*x1ry1r5.*y1r1.*y2r2-x2r5.*x1ry1r1.*y1r3.*y2r2-x2r5.*x1ry1r2.*y1r1.*y2r3-x2r5.*x1ry1r3.*y1r2.*y2r1-x2r1.*x1ry1r2.*y1r5.*y2r4-x2r1.*x1ry1r4.*y1r2.*y2r5-x2r1.*x1ry1r5.*y1r4.*y2r2-x2r2.*x1ry1r1.*y1r4.*y2r5-x2r2.*x1ry1r4.*y1r5.*y2r1-x2r2.*x1ry1r5.*y1r1.*y2r4-x2r4.*x1ry1r1.*y1r5.*y2r2-x2r4.*x1ry1r2.*y1r1.*y2r5-x2r4.*x1ry1r5.*y1r2.*y2r1-x2r5.*x1ry1r1.*y1r2.*y2r4-x2r5.*x1ry1r2.*y1r4.*y2r1-x2r5.*x1ry1r4.*y1r1.*y2r2-x2r1.*x1ry1r3.*y1r4.*y2r5-x2r1.*x1ry1r4.*y1r5.*y2r3-x2r1.*x1ry1r5.*y1r3.*y2r4-x2r3.*x1ry1r1.*y1r5.*y2r4-x2r3.*x1ry1r4.*y1r1.*y2r5-x2r3.*x1ry1r5.*y1r4.*y2r1-x2r4.*x1ry1r1.*y1r3.*y2r5-x2r4.*x1ry1r3.*y1r5.*y2r1-x2r4.*x1ry1r5.*y1r1.*y2r3-x2r5.*x1ry1r1.*y1r4.*y2r3-x2r5.*x1ry1r3.*y1r1.*y2r4-x2r5.*x1ry1r4.*y1r3.*y2r1-x2r2.*x1ry1r3.*y1r5.*y2r4-x2r2.*x1ry1r4.*y1r3.*y2r5-x2r2.*x1ry1r5.*y1r4.*y2r3-x2r3.*x1ry1r2.*y1r4.*y2r5-x2r3.*x1ry1r4.*y1r5.*y2r2-x2r3.*x1ry1r5.*y1r2.*y2r4-x2r4.*x1ry1r2.*y1r5.*y2r3-x2r4.*x1ry1r3.*y1r2.*y2r5-x2r4.*x1ry1r5.*y1r3.*y2r2-x2r5.*x1ry1r2.*y1r3.*y2r4-x2r5.*x1ry1r3.*y1r4.*y2r2-x2r5.*x1ry1r4.*y1r2.*y2r3),t75.*(t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13-x2r2.*x1ry1r3.*y1r5.*y2r4-x2r2.*x1ry1r4.*y1r3.*y2r5-x2r2.*x1ry1r5.*y1r4.*y2r3-x2r3.*x1ry1r2.*y1r4.*y2r5-x2r3.*x1ry1r4.*y1r5.*y2r2-x2r3.*x1ry1r5.*y1r2.*y2r4-x2r4.*x1ry1r2.*y1r5.*y2r3-x2r4.*x1ry1r3.*y1r2.*y2r5-x2r4.*x1ry1r5.*y1r3.*y2r2-x2r5.*x1ry1r2.*y1r3.*y2r4-x2r5.*x1ry1r3.*y1r4.*y2r2-x2r5.*x1ry1r4.*y1r2.*y2r3),t75.*(t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t86+t87-x2r1.*x1ry1r3.*y1r4.*y2r5-x2r1.*x1ry1r4.*y1r5.*y2r3-x2r1.*x1ry1r5.*y1r3.*y2r4-x2r3.*x1ry1r1.*y1r5.*y2r4-x2r3.*x1ry1r4.*y1r1.*y2r5-x2r3.*x1ry1r5.*y1r4.*y2r1-x2r4.*x1ry1r1.*y1r3.*y2r5-x2r4.*x1ry1r3.*y1r5.*y2r1-x2r4.*x1ry1r5.*y1r1.*y2r3-x2r5.*x1ry1r1.*y1r4.*y2r3-x2r5.*x1ry1r3.*y1r1.*y2r4-x2r5.*x1ry1r4.*y1r3.*y2r1),t75.*(t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159-x2r1.*x1ry1r2.*y1r5.*y2r4-x2r1.*x1ry1r4.*y1r2.*y2r5-x2r1.*x1ry1r5.*y1r4.*y2r2-x2r2.*x1ry1r1.*y1r4.*y2r5-x2r2.*x1ry1r4.*y1r5.*y2r1-x2r2.*x1ry1r5.*y1r1.*y2r4-x2r4.*x1ry1r1.*y1r5.*y2r2-x2r4.*x1ry1r2.*y1r1.*y2r5-x2r4.*x1ry1r5.*y1r2.*y2r1-x2r5.*x1ry1r1.*y1r2.*y2r4-x2r5.*x1ry1r2.*y1r4.*y2r1-x2r5.*x1ry1r4.*y1r1.*y2r2),t75.*(t160+t161+t162+t163+t164+t165+t166+t167+t168+t169+t170+t171-x2r1.*x1ry1r2.*y1r3.*y2r5-x2r1.*x1ry1r3.*y1r5.*y2r2-x2r1.*x1ry1r5.*y1r2.*y2r3-x2r2.*x1ry1r1.*y1r5.*y2r3-x2r2.*x1ry1r3.*y1r1.*y2r5-x2r2.*x1ry1r5.*y1r3.*y2r1-x2r3.*x1ry1r1.*y1r2.*y2r5-x2r3.*x1ry1r2.*y1r5.*y2r1-x2r3.*x1ry1r5.*y1r1.*y2r2-x2r5.*x1ry1r1.*y1r3.*y2r2-x2r5.*x1ry1r2.*y1r1.*y2r3-x2r5.*x1ry1r3.*y1r2.*y2r1),t75.*(t172+t173+t174+t175+t176+t177+t178+t179+t180+t181+t182+t183-x2r1.*x1ry1r2.*y1r4.*y2r3-x2r1.*x1ry1r3.*y1r2.*y2r4-x2r1.*x1ry1r4.*y1r3.*y2r2-x2r2.*x1ry1r1.*y1r3.*y2r4-x2r2.*x1ry1r3.*y1r4.*y2r1-x2r2.*x1ry1r4.*y1r1.*y2r3-x2r3.*x1ry1r1.*y1r4.*y2r2-x2r3.*x1ry1r2.*y1r1.*y2r4-x2r3.*x1ry1r4.*y1r2.*y2r1-x2r4.*x1ry1r1.*y1r2.*y2r3-x2r4.*x1ry1r2.*y1r3.*y2r1-x2r4.*x1ry1r3.*y1r1.*y2r2)];
if nargout > 1
    t184 = x1r2.*x2r3.*x1ry1r4.*y2r5;
    t185 = x1r2.*x2r4.*x1ry1r5.*y2r3;
    t186 = x1r2.*x2r5.*x1ry1r3.*y2r4;
    t187 = x1r3.*x2r2.*x1ry1r5.*y2r4;
    t188 = x1r3.*x2r4.*x1ry1r2.*y2r5;
    t189 = x1r3.*x2r5.*x1ry1r4.*y2r2;
    t190 = x1r4.*x2r2.*x1ry1r3.*y2r5;
    t191 = x1r4.*x2r3.*x1ry1r5.*y2r2;
    t192 = x1r4.*x2r5.*x1ry1r2.*y2r3;
    t193 = x1r5.*x2r2.*x1ry1r4.*y2r3;
    t194 = x1r5.*x2r3.*x1ry1r2.*y2r4;
    t195 = x1r5.*x2r4.*x1ry1r3.*y2r2;
    t196 = x1r1.*x2r3.*x1ry1r5.*y2r4;
    t197 = x1r1.*x2r4.*x1ry1r3.*y2r5;
    t198 = x1r1.*x2r5.*x1ry1r4.*y2r3;
    t199 = x1r3.*x2r1.*x1ry1r4.*y2r5;
    t200 = x1r3.*x2r4.*x1ry1r5.*y2r1;
    t201 = x1r3.*x2r5.*x1ry1r1.*y2r4;
    t202 = x1r4.*x2r1.*x1ry1r5.*y2r3;
    t203 = x1r4.*x2r3.*x1ry1r1.*y2r5;
    t204 = x1r4.*x2r5.*x1ry1r3.*y2r1;
    t205 = x1r5.*x2r1.*x1ry1r3.*y2r4;
    t206 = x1r5.*x2r3.*x1ry1r4.*y2r1;
    t207 = x1r5.*x2r4.*x1ry1r1.*y2r3;
    t208 = x1r1.*x2r2.*x1ry1r4.*y2r5;
    t209 = x1r1.*x2r4.*x1ry1r5.*y2r2;
    t210 = x1r1.*x2r5.*x1ry1r2.*y2r4;
    t211 = x1r2.*x2r1.*x1ry1r5.*y2r4;
    t212 = x1r2.*x2r4.*x1ry1r1.*y2r5;
    t213 = x1r2.*x2r5.*x1ry1r4.*y2r1;
    t214 = x1r4.*x2r1.*x1ry1r2.*y2r5;
    t215 = x1r4.*x2r2.*x1ry1r5.*y2r1;
    t216 = x1r4.*x2r5.*x1ry1r1.*y2r2;
    t217 = x1r5.*x2r1.*x1ry1r4.*y2r2;
    t218 = x1r5.*x2r2.*x1ry1r1.*y2r4;
    t219 = x1r5.*x2r4.*x1ry1r2.*y2r1;
    t220 = x1r1.*x2r2.*x1ry1r5.*y2r3;
    t221 = x1r1.*x2r3.*x1ry1r2.*y2r5;
    t222 = x1r1.*x2r5.*x1ry1r3.*y2r2;
    t223 = x1r2.*x2r1.*x1ry1r3.*y2r5;
    t224 = x1r2.*x2r3.*x1ry1r5.*y2r1;
    t225 = x1r2.*x2r5.*x1ry1r1.*y2r3;
    t226 = x1r3.*x2r1.*x1ry1r5.*y2r2;
    t227 = x1r3.*x2r2.*x1ry1r1.*y2r5;
    t228 = x1r3.*x2r5.*x1ry1r2.*y2r1;
    t229 = x1r5.*x2r1.*x1ry1r2.*y2r3;
    t230 = x1r5.*x2r2.*x1ry1r3.*y2r1;
    t231 = x1r5.*x2r3.*x1ry1r1.*y2r2;
    t232 = x1r1.*x2r2.*x1ry1r3.*y2r4;
    t233 = x1r1.*x2r3.*x1ry1r4.*y2r2;
    t234 = x1r1.*x2r4.*x1ry1r2.*y2r3;
    t235 = x1r2.*x2r1.*x1ry1r4.*y2r3;
    t236 = x1r2.*x2r3.*x1ry1r1.*y2r4;
    t237 = x1r2.*x2r4.*x1ry1r3.*y2r1;
    t238 = x1r3.*x2r1.*x1ry1r2.*y2r4;
    t239 = x1r3.*x2r2.*x1ry1r4.*y2r1;
    t240 = x1r3.*x2r4.*x1ry1r1.*y2r2;
    t241 = x1r4.*x2r1.*x1ry1r3.*y2r2;
    t242 = x1r4.*x2r2.*x1ry1r1.*y2r3;
    t243 = x1r4.*x2r3.*x1ry1r2.*y2r1;
    duyva = [t75.*(t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t200+t201+t202+t203+t204+t205+t206+t207+t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+t218+t219+t220+t221+t222+t223+t224+t225+t226+t227+t228+t229+t230+t231+t232+t233+t234+t235+t236+t237+t238+t239+t240+t241+t242+t243-x1r1.*x2r2.*x1ry1r4.*y2r3-x1r1.*x2r3.*x1ry1r2.*y2r4-x1r1.*x2r4.*x1ry1r3.*y2r2-x1r2.*x2r1.*x1ry1r3.*y2r4-x1r2.*x2r3.*x1ry1r4.*y2r1-x1r2.*x2r4.*x1ry1r1.*y2r3-x1r3.*x2r1.*x1ry1r4.*y2r2-x1r3.*x2r2.*x1ry1r1.*y2r4-x1r3.*x2r4.*x1ry1r2.*y2r1-x1r4.*x2r1.*x1ry1r2.*y2r3-x1r4.*x2r2.*x1ry1r3.*y2r1-x1r4.*x2r3.*x1ry1r1.*y2r2-x1r1.*x2r2.*x1ry1r3.*y2r5-x1r1.*x2r3.*x1ry1r5.*y2r2-x1r1.*x2r5.*x1ry1r2.*y2r3-x1r2.*x2r1.*x1ry1r5.*y2r3-x1r2.*x2r3.*x1ry1r1.*y2r5-x1r2.*x2r5.*x1ry1r3.*y2r1-x1r3.*x2r1.*x1ry1r2.*y2r5-x1r3.*x2r2.*x1ry1r5.*y2r1-x1r3.*x2r5.*x1ry1r1.*y2r2-x1r5.*x2r1.*x1ry1r3.*y2r2-x1r5.*x2r2.*x1ry1r1.*y2r3-x1r5.*x2r3.*x1ry1r2.*y2r1-x1r1.*x2r2.*x1ry1r5.*y2r4-x1r1.*x2r4.*x1ry1r2.*y2r5-x1r1.*x2r5.*x1ry1r4.*y2r2-x1r2.*x2r1.*x1ry1r4.*y2r5-x1r2.*x2r4.*x1ry1r5.*y2r1-x1r2.*x2r5.*x1ry1r1.*y2r4-x1r4.*x2r1.*x1ry1r5.*y2r2-x1r4.*x2r2.*x1ry1r1.*y2r5-x1r4.*x2r5.*x1ry1r2.*y2r1-x1r5.*x2r1.*x1ry1r2.*y2r4-x1r5.*x2r2.*x1ry1r4.*y2r1-x1r5.*x2r4.*x1ry1r1.*y2r2-x1r1.*x2r3.*x1ry1r4.*y2r5-x1r1.*x2r4.*x1ry1r5.*y2r3-x1r1.*x2r5.*x1ry1r3.*y2r4-x1r3.*x2r1.*x1ry1r5.*y2r4-x1r3.*x2r4.*x1ry1r1.*y2r5-x1r3.*x2r5.*x1ry1r4.*y2r1-x1r4.*x2r1.*x1ry1r3.*y2r5-x1r4.*x2r3.*x1ry1r5.*y2r1-x1r4.*x2r5.*x1ry1r1.*y2r3-x1r5.*x2r1.*x1ry1r4.*y2r3-x1r5.*x2r3.*x1ry1r1.*y2r4-x1r5.*x2r4.*x1ry1r3.*y2r1-x1r2.*x2r3.*x1ry1r5.*y2r4-x1r2.*x2r4.*x1ry1r3.*y2r5-x1r2.*x2r5.*x1ry1r4.*y2r3-x1r3.*x2r2.*x1ry1r4.*y2r5-x1r3.*x2r4.*x1ry1r5.*y2r2-x1r3.*x2r5.*x1ry1r2.*y2r4-x1r4.*x2r2.*x1ry1r5.*y2r3-x1r4.*x2r3.*x1ry1r2.*y2r5-x1r4.*x2r5.*x1ry1r3.*y2r2-x1r5.*x2r2.*x1ry1r3.*y2r4-x1r5.*x2r3.*x1ry1r4.*y2r2-x1r5.*x2r4.*x1ry1r2.*y2r3),-t75.*(t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t194+t195-x1r2.*x2r3.*x1ry1r5.*y2r4-x1r2.*x2r4.*x1ry1r3.*y2r5-x1r2.*x2r5.*x1ry1r4.*y2r3-x1r3.*x2r2.*x1ry1r4.*y2r5-x1r3.*x2r4.*x1ry1r5.*y2r2-x1r3.*x2r5.*x1ry1r2.*y2r4-x1r4.*x2r2.*x1ry1r5.*y2r3-x1r4.*x2r3.*x1ry1r2.*y2r5-x1r4.*x2r5.*x1ry1r3.*y2r2-x1r5.*x2r2.*x1ry1r3.*y2r4-x1r5.*x2r3.*x1ry1r4.*y2r2-x1r5.*x2r4.*x1ry1r2.*y2r3),-t75.*(t196+t197+t198+t199+t200+t201+t202+t203+t204+t205+t206+t207-x1r1.*x2r3.*x1ry1r4.*y2r5-x1r1.*x2r4.*x1ry1r5.*y2r3-x1r1.*x2r5.*x1ry1r3.*y2r4-x1r3.*x2r1.*x1ry1r5.*y2r4-x1r3.*x2r4.*x1ry1r1.*y2r5-x1r3.*x2r5.*x1ry1r4.*y2r1-x1r4.*x2r1.*x1ry1r3.*y2r5-x1r4.*x2r3.*x1ry1r5.*y2r1-x1r4.*x2r5.*x1ry1r1.*y2r3-x1r5.*x2r1.*x1ry1r4.*y2r3-x1r5.*x2r3.*x1ry1r1.*y2r4-x1r5.*x2r4.*x1ry1r3.*y2r1),-t75.*(t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+t218+t219-x1r1.*x2r2.*x1ry1r5.*y2r4-x1r1.*x2r4.*x1ry1r2.*y2r5-x1r1.*x2r5.*x1ry1r4.*y2r2-x1r2.*x2r1.*x1ry1r4.*y2r5-x1r2.*x2r4.*x1ry1r5.*y2r1-x1r2.*x2r5.*x1ry1r1.*y2r4-x1r4.*x2r1.*x1ry1r5.*y2r2-x1r4.*x2r2.*x1ry1r1.*y2r5-x1r4.*x2r5.*x1ry1r2.*y2r1-x1r5.*x2r1.*x1ry1r2.*y2r4-x1r5.*x2r2.*x1ry1r4.*y2r1-x1r5.*x2r4.*x1ry1r1.*y2r2),-t75.*(t220+t221+t222+t223+t224+t225+t226+t227+t228+t229+t230+t231-x1r1.*x2r2.*x1ry1r3.*y2r5-x1r1.*x2r3.*x1ry1r5.*y2r2-x1r1.*x2r5.*x1ry1r2.*y2r3-x1r2.*x2r1.*x1ry1r5.*y2r3-x1r2.*x2r3.*x1ry1r1.*y2r5-x1r2.*x2r5.*x1ry1r3.*y2r1-x1r3.*x2r1.*x1ry1r2.*y2r5-x1r3.*x2r2.*x1ry1r5.*y2r1-x1r3.*x2r5.*x1ry1r1.*y2r2-x1r5.*x2r1.*x1ry1r3.*y2r2-x1r5.*x2r2.*x1ry1r1.*y2r3-x1r5.*x2r3.*x1ry1r2.*y2r1),-t75.*(t232+t233+t234+t235+t236+t237+t238+t239+t240+t241+t242+t243-x1r1.*x2r2.*x1ry1r4.*y2r3-x1r1.*x2r3.*x1ry1r2.*y2r4-x1r1.*x2r4.*x1ry1r3.*y2r2-x1r2.*x2r1.*x1ry1r3.*y2r4-x1r2.*x2r3.*x1ry1r4.*y2r1-x1r2.*x2r4.*x1ry1r1.*y2r3-x1r3.*x2r1.*x1ry1r4.*y2r2-x1r3.*x2r2.*x1ry1r1.*y2r4-x1r3.*x2r4.*x1ry1r2.*y2r1-x1r4.*x2r1.*x1ry1r2.*y2r3-x1r4.*x2r2.*x1ry1r3.*y2r1-x1r4.*x2r3.*x1ry1r1.*y2r2)];
end
if nargout > 2
    t244 = x1r2.*x1ry1r3.*y1r4.*y2r5;
    t245 = x1r2.*x1ry1r4.*y1r5.*y2r3;
    t246 = x1r2.*x1ry1r5.*y1r3.*y2r4;
    t247 = x1r3.*x1ry1r2.*y1r5.*y2r4;
    t248 = x1r3.*x1ry1r4.*y1r2.*y2r5;
    t249 = x1r3.*x1ry1r5.*y1r4.*y2r2;
    t250 = x1r4.*x1ry1r2.*y1r3.*y2r5;
    t251 = x1r4.*x1ry1r3.*y1r5.*y2r2;
    t252 = x1r4.*x1ry1r5.*y1r2.*y2r3;
    t253 = x1r5.*x1ry1r2.*y1r4.*y2r3;
    t254 = x1r5.*x1ry1r3.*y1r2.*y2r4;
    t255 = x1r5.*x1ry1r4.*y1r3.*y2r2;
    t256 = x1r1.*x1ry1r3.*y1r5.*y2r4;
    t257 = x1r1.*x1ry1r4.*y1r3.*y2r5;
    t258 = x1r1.*x1ry1r5.*y1r4.*y2r3;
    t259 = x1r3.*x1ry1r1.*y1r4.*y2r5;
    t260 = x1r3.*x1ry1r4.*y1r5.*y2r1;
    t261 = x1r3.*x1ry1r5.*y1r1.*y2r4;
    t262 = x1r4.*x1ry1r1.*y1r5.*y2r3;
    t263 = x1r4.*x1ry1r3.*y1r1.*y2r5;
    t264 = x1r4.*x1ry1r5.*y1r3.*y2r1;
    t265 = x1r5.*x1ry1r1.*y1r3.*y2r4;
    t266 = x1r5.*x1ry1r3.*y1r4.*y2r1;
    t267 = x1r5.*x1ry1r4.*y1r1.*y2r3;
    t268 = x1r1.*x1ry1r2.*y1r4.*y2r5;
    t269 = x1r1.*x1ry1r4.*y1r5.*y2r2;
    t270 = x1r1.*x1ry1r5.*y1r2.*y2r4;
    t271 = x1r2.*x1ry1r1.*y1r5.*y2r4;
    t272 = x1r2.*x1ry1r4.*y1r1.*y2r5;
    t273 = x1r2.*x1ry1r5.*y1r4.*y2r1;
    t274 = x1r4.*x1ry1r1.*y1r2.*y2r5;
    t275 = x1r4.*x1ry1r2.*y1r5.*y2r1;
    t276 = x1r4.*x1ry1r5.*y1r1.*y2r2;
    t277 = x1r5.*x1ry1r1.*y1r4.*y2r2;
    t278 = x1r5.*x1ry1r2.*y1r1.*y2r4;
    t279 = x1r5.*x1ry1r4.*y1r2.*y2r1;
    t280 = x1r1.*x1ry1r2.*y1r5.*y2r3;
    t281 = x1r1.*x1ry1r3.*y1r2.*y2r5;
    t282 = x1r1.*x1ry1r5.*y1r3.*y2r2;
    t283 = x1r2.*x1ry1r1.*y1r3.*y2r5;
    t284 = x1r2.*x1ry1r3.*y1r5.*y2r1;
    t285 = x1r2.*x1ry1r5.*y1r1.*y2r3;
    t286 = x1r3.*x1ry1r1.*y1r5.*y2r2;
    t287 = x1r3.*x1ry1r2.*y1r1.*y2r5;
    t288 = x1r3.*x1ry1r5.*y1r2.*y2r1;
    t289 = x1r5.*x1ry1r1.*y1r2.*y2r3;
    t290 = x1r5.*x1ry1r2.*y1r3.*y2r1;
    t291 = x1r5.*x1ry1r3.*y1r1.*y2r2;
    t292 = x1r1.*x1ry1r2.*y1r3.*y2r4;
    t293 = x1r1.*x1ry1r3.*y1r4.*y2r2;
    t294 = x1r1.*x1ry1r4.*y1r2.*y2r3;
    t295 = x1r2.*x1ry1r1.*y1r4.*y2r3;
    t296 = x1r2.*x1ry1r3.*y1r1.*y2r4;
    t297 = x1r2.*x1ry1r4.*y1r3.*y2r1;
    t298 = x1r3.*x1ry1r1.*y1r2.*y2r4;
    t299 = x1r3.*x1ry1r2.*y1r4.*y2r1;
    t300 = x1r3.*x1ry1r4.*y1r1.*y2r2;
    t301 = x1r4.*x1ry1r1.*y1r3.*y2r2;
    t302 = x1r4.*x1ry1r2.*y1r1.*y2r3;
    t303 = x1r4.*x1ry1r3.*y1r2.*y2r1;
    duxxva = [t75.*(t244+t245+t246+t247+t248+t249+t250+t251+t252+t253+t254+t255+t256+t257+t258+t259+t260+t261+t262+t263+t264+t265+t266+t267+t268+t269+t270+t271+t272+t273+t274+t275+t276+t277+t278+t279+t280+t281+t282+t283+t284+t285+t286+t287+t288+t289+t290+t291+t292+t293+t294+t295+t296+t297+t298+t299+t300+t301+t302+t303-x1r1.*x1ry1r2.*y1r4.*y2r3-x1r1.*x1ry1r3.*y1r2.*y2r4-x1r1.*x1ry1r4.*y1r3.*y2r2-x1r2.*x1ry1r1.*y1r3.*y2r4-x1r2.*x1ry1r3.*y1r4.*y2r1-x1r2.*x1ry1r4.*y1r1.*y2r3-x1r3.*x1ry1r1.*y1r4.*y2r2-x1r3.*x1ry1r2.*y1r1.*y2r4-x1r3.*x1ry1r4.*y1r2.*y2r1-x1r4.*x1ry1r1.*y1r2.*y2r3-x1r4.*x1ry1r2.*y1r3.*y2r1-x1r4.*x1ry1r3.*y1r1.*y2r2-x1r1.*x1ry1r2.*y1r3.*y2r5-x1r1.*x1ry1r3.*y1r5.*y2r2-x1r1.*x1ry1r5.*y1r2.*y2r3-x1r2.*x1ry1r1.*y1r5.*y2r3-x1r2.*x1ry1r3.*y1r1.*y2r5-x1r2.*x1ry1r5.*y1r3.*y2r1-x1r3.*x1ry1r1.*y1r2.*y2r5-x1r3.*x1ry1r2.*y1r5.*y2r1-x1r3.*x1ry1r5.*y1r1.*y2r2-x1r5.*x1ry1r1.*y1r3.*y2r2-x1r5.*x1ry1r2.*y1r1.*y2r3-x1r5.*x1ry1r3.*y1r2.*y2r1-x1r1.*x1ry1r2.*y1r5.*y2r4-x1r1.*x1ry1r4.*y1r2.*y2r5-x1r1.*x1ry1r5.*y1r4.*y2r2-x1r2.*x1ry1r1.*y1r4.*y2r5-x1r2.*x1ry1r4.*y1r5.*y2r1-x1r2.*x1ry1r5.*y1r1.*y2r4-x1r4.*x1ry1r1.*y1r5.*y2r2-x1r4.*x1ry1r2.*y1r1.*y2r5-x1r4.*x1ry1r5.*y1r2.*y2r1-x1r5.*x1ry1r1.*y1r2.*y2r4-x1r5.*x1ry1r2.*y1r4.*y2r1-x1r5.*x1ry1r4.*y1r1.*y2r2-x1r1.*x1ry1r3.*y1r4.*y2r5-x1r1.*x1ry1r4.*y1r5.*y2r3-x1r1.*x1ry1r5.*y1r3.*y2r4-x1r3.*x1ry1r1.*y1r5.*y2r4-x1r3.*x1ry1r4.*y1r1.*y2r5-x1r3.*x1ry1r5.*y1r4.*y2r1-x1r4.*x1ry1r1.*y1r3.*y2r5-x1r4.*x1ry1r3.*y1r5.*y2r1-x1r4.*x1ry1r5.*y1r1.*y2r3-x1r5.*x1ry1r1.*y1r4.*y2r3-x1r5.*x1ry1r3.*y1r1.*y2r4-x1r5.*x1ry1r4.*y1r3.*y2r1-x1r2.*x1ry1r3.*y1r5.*y2r4-x1r2.*x1ry1r4.*y1r3.*y2r5-x1r2.*x1ry1r5.*y1r4.*y2r3-x1r3.*x1ry1r2.*y1r4.*y2r5-x1r3.*x1ry1r4.*y1r5.*y2r2-x1r3.*x1ry1r5.*y1r2.*y2r4-x1r4.*x1ry1r2.*y1r5.*y2r3-x1r4.*x1ry1r3.*y1r2.*y2r5-x1r4.*x1ry1r5.*y1r3.*y2r2-x1r5.*x1ry1r2.*y1r3.*y2r4-x1r5.*x1ry1r3.*y1r4.*y2r2-x1r5.*x1ry1r4.*y1r2.*y2r3).*2.0,t75.*(t244+t245+t246+t247+t248+t249+t250+t251+t252+t253+t254+t255-x1r2.*x1ry1r3.*y1r5.*y2r4-x1r2.*x1ry1r4.*y1r3.*y2r5-x1r2.*x1ry1r5.*y1r4.*y2r3-x1r3.*x1ry1r2.*y1r4.*y2r5-x1r3.*x1ry1r4.*y1r5.*y2r2-x1r3.*x1ry1r5.*y1r2.*y2r4-x1r4.*x1ry1r2.*y1r5.*y2r3-x1r4.*x1ry1r3.*y1r2.*y2r5-x1r4.*x1ry1r5.*y1r3.*y2r2-x1r5.*x1ry1r2.*y1r3.*y2r4-x1r5.*x1ry1r3.*y1r4.*y2r2-x1r5.*x1ry1r4.*y1r2.*y2r3).*-2.0,t75.*(t256+t257+t258+t259+t260+t261+t262+t263+t264+t265+t266+t267-x1r1.*x1ry1r3.*y1r4.*y2r5-x1r1.*x1ry1r4.*y1r5.*y2r3-x1r1.*x1ry1r5.*y1r3.*y2r4-x1r3.*x1ry1r1.*y1r5.*y2r4-x1r3.*x1ry1r4.*y1r1.*y2r5-x1r3.*x1ry1r5.*y1r4.*y2r1-x1r4.*x1ry1r1.*y1r3.*y2r5-x1r4.*x1ry1r3.*y1r5.*y2r1-x1r4.*x1ry1r5.*y1r1.*y2r3-x1r5.*x1ry1r1.*y1r4.*y2r3-x1r5.*x1ry1r3.*y1r1.*y2r4-x1r5.*x1ry1r4.*y1r3.*y2r1).*-2.0,t75.*(t268+t269+t270+t271+t272+t273+t274+t275+t276+t277+t278+t279-x1r1.*x1ry1r2.*y1r5.*y2r4-x1r1.*x1ry1r4.*y1r2.*y2r5-x1r1.*x1ry1r5.*y1r4.*y2r2-x1r2.*x1ry1r1.*y1r4.*y2r5-x1r2.*x1ry1r4.*y1r5.*y2r1-x1r2.*x1ry1r5.*y1r1.*y2r4-x1r4.*x1ry1r1.*y1r5.*y2r2-x1r4.*x1ry1r2.*y1r1.*y2r5-x1r4.*x1ry1r5.*y1r2.*y2r1-x1r5.*x1ry1r1.*y1r2.*y2r4-x1r5.*x1ry1r2.*y1r4.*y2r1-x1r5.*x1ry1r4.*y1r1.*y2r2).*-2.0,t75.*(t280+t281+t282+t283+t284+t285+t286+t287+t288+t289+t290+t291-x1r1.*x1ry1r2.*y1r3.*y2r5-x1r1.*x1ry1r3.*y1r5.*y2r2-x1r1.*x1ry1r5.*y1r2.*y2r3-x1r2.*x1ry1r1.*y1r5.*y2r3-x1r2.*x1ry1r3.*y1r1.*y2r5-x1r2.*x1ry1r5.*y1r3.*y2r1-x1r3.*x1ry1r1.*y1r2.*y2r5-x1r3.*x1ry1r2.*y1r5.*y2r1-x1r3.*x1ry1r5.*y1r1.*y2r2-x1r5.*x1ry1r1.*y1r3.*y2r2-x1r5.*x1ry1r2.*y1r1.*y2r3-x1r5.*x1ry1r3.*y1r2.*y2r1).*-2.0,t75.*(t292+t293+t294+t295+t296+t297+t298+t299+t300+t301+t302+t303-x1r1.*x1ry1r2.*y1r4.*y2r3-x1r1.*x1ry1r3.*y1r2.*y2r4-x1r1.*x1ry1r4.*y1r3.*y2r2-x1r2.*x1ry1r1.*y1r3.*y2r4-x1r2.*x1ry1r3.*y1r4.*y2r1-x1r2.*x1ry1r4.*y1r1.*y2r3-x1r3.*x1ry1r1.*y1r4.*y2r2-x1r3.*x1ry1r2.*y1r1.*y2r4-x1r3.*x1ry1r4.*y1r2.*y2r1-x1r4.*x1ry1r1.*y1r2.*y2r3-x1r4.*x1ry1r2.*y1r3.*y2r1-x1r4.*x1ry1r3.*y1r1.*y2r2).*-2.0];
end
if nargout > 3
    t304 = x1r2.*x2r3.*x1ry1r4.*y1r5;
    t305 = x1r2.*x2r4.*x1ry1r5.*y1r3;
    t306 = x1r2.*x2r5.*x1ry1r3.*y1r4;
    t307 = x1r3.*x2r2.*x1ry1r5.*y1r4;
    t308 = x1r3.*x2r4.*x1ry1r2.*y1r5;
    t309 = x1r3.*x2r5.*x1ry1r4.*y1r2;
    t310 = x1r4.*x2r2.*x1ry1r3.*y1r5;
    t311 = x1r4.*x2r3.*x1ry1r5.*y1r2;
    t312 = x1r4.*x2r5.*x1ry1r2.*y1r3;
    t313 = x1r5.*x2r2.*x1ry1r4.*y1r3;
    t314 = x1r5.*x2r3.*x1ry1r2.*y1r4;
    t315 = x1r5.*x2r4.*x1ry1r3.*y1r2;
    t316 = x1r1.*x2r3.*x1ry1r5.*y1r4;
    t317 = x1r1.*x2r4.*x1ry1r3.*y1r5;
    t318 = x1r1.*x2r5.*x1ry1r4.*y1r3;
    t319 = x1r3.*x2r1.*x1ry1r4.*y1r5;
    t320 = x1r3.*x2r4.*x1ry1r5.*y1r1;
    t321 = x1r3.*x2r5.*x1ry1r1.*y1r4;
    t322 = x1r4.*x2r1.*x1ry1r5.*y1r3;
    t323 = x1r4.*x2r3.*x1ry1r1.*y1r5;
    t324 = x1r4.*x2r5.*x1ry1r3.*y1r1;
    t325 = x1r5.*x2r1.*x1ry1r3.*y1r4;
    t326 = x1r5.*x2r3.*x1ry1r4.*y1r1;
    t327 = x1r5.*x2r4.*x1ry1r1.*y1r3;
    t328 = x1r1.*x2r2.*x1ry1r4.*y1r5;
    t329 = x1r1.*x2r4.*x1ry1r5.*y1r2;
    t330 = x1r1.*x2r5.*x1ry1r2.*y1r4;
    t331 = x1r2.*x2r1.*x1ry1r5.*y1r4;
    t332 = x1r2.*x2r4.*x1ry1r1.*y1r5;
    t333 = x1r2.*x2r5.*x1ry1r4.*y1r1;
    t334 = x1r4.*x2r1.*x1ry1r2.*y1r5;
    t335 = x1r4.*x2r2.*x1ry1r5.*y1r1;
    t336 = x1r4.*x2r5.*x1ry1r1.*y1r2;
    t337 = x1r5.*x2r1.*x1ry1r4.*y1r2;
    t338 = x1r5.*x2r2.*x1ry1r1.*y1r4;
    t339 = x1r5.*x2r4.*x1ry1r2.*y1r1;
    t340 = x1r1.*x2r2.*x1ry1r5.*y1r3;
    t341 = x1r1.*x2r3.*x1ry1r2.*y1r5;
    t342 = x1r1.*x2r5.*x1ry1r3.*y1r2;
    t343 = x1r2.*x2r1.*x1ry1r3.*y1r5;
    t344 = x1r2.*x2r3.*x1ry1r5.*y1r1;
    t345 = x1r2.*x2r5.*x1ry1r1.*y1r3;
    t346 = x1r3.*x2r1.*x1ry1r5.*y1r2;
    t347 = x1r3.*x2r2.*x1ry1r1.*y1r5;
    t348 = x1r3.*x2r5.*x1ry1r2.*y1r1;
    t349 = x1r5.*x2r1.*x1ry1r2.*y1r3;
    t350 = x1r5.*x2r2.*x1ry1r3.*y1r1;
    t351 = x1r5.*x2r3.*x1ry1r1.*y1r2;
    t352 = x1r1.*x2r2.*x1ry1r3.*y1r4;
    t353 = x1r1.*x2r3.*x1ry1r4.*y1r2;
    t354 = x1r1.*x2r4.*x1ry1r2.*y1r3;
    t355 = x1r2.*x2r1.*x1ry1r4.*y1r3;
    t356 = x1r2.*x2r3.*x1ry1r1.*y1r4;
    t357 = x1r2.*x2r4.*x1ry1r3.*y1r1;
    t358 = x1r3.*x2r1.*x1ry1r2.*y1r4;
    t359 = x1r3.*x2r2.*x1ry1r4.*y1r1;
    t360 = x1r3.*x2r4.*x1ry1r1.*y1r2;
    t361 = x1r4.*x2r1.*x1ry1r3.*y1r2;
    t362 = x1r4.*x2r2.*x1ry1r1.*y1r3;
    t363 = x1r4.*x2r3.*x1ry1r2.*y1r1;
    duyyva = [t75.*(t304+t305+t306+t307+t308+t309+t310+t311+t312+t313+t314+t315+t316+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t331+t332+t333+t334+t335+t336+t337+t338+t339+t340+t341+t342+t343+t344+t345+t346+t347+t348+t349+t350+t351+t352+t353+t354+t355+t356+t357+t358+t359+t360+t361+t362+t363-x1r1.*x2r2.*x1ry1r4.*y1r3-x1r1.*x2r3.*x1ry1r2.*y1r4-x1r1.*x2r4.*x1ry1r3.*y1r2-x1r2.*x2r1.*x1ry1r3.*y1r4-x1r2.*x2r3.*x1ry1r4.*y1r1-x1r2.*x2r4.*x1ry1r1.*y1r3-x1r3.*x2r1.*x1ry1r4.*y1r2-x1r3.*x2r2.*x1ry1r1.*y1r4-x1r3.*x2r4.*x1ry1r2.*y1r1-x1r4.*x2r1.*x1ry1r2.*y1r3-x1r4.*x2r2.*x1ry1r3.*y1r1-x1r4.*x2r3.*x1ry1r1.*y1r2-x1r1.*x2r2.*x1ry1r3.*y1r5-x1r1.*x2r3.*x1ry1r5.*y1r2-x1r1.*x2r5.*x1ry1r2.*y1r3-x1r2.*x2r1.*x1ry1r5.*y1r3-x1r2.*x2r3.*x1ry1r1.*y1r5-x1r2.*x2r5.*x1ry1r3.*y1r1-x1r3.*x2r1.*x1ry1r2.*y1r5-x1r3.*x2r2.*x1ry1r5.*y1r1-x1r3.*x2r5.*x1ry1r1.*y1r2-x1r5.*x2r1.*x1ry1r3.*y1r2-x1r5.*x2r2.*x1ry1r1.*y1r3-x1r5.*x2r3.*x1ry1r2.*y1r1-x1r1.*x2r2.*x1ry1r5.*y1r4-x1r1.*x2r4.*x1ry1r2.*y1r5-x1r1.*x2r5.*x1ry1r4.*y1r2-x1r2.*x2r1.*x1ry1r4.*y1r5-x1r2.*x2r4.*x1ry1r5.*y1r1-x1r2.*x2r5.*x1ry1r1.*y1r4-x1r4.*x2r1.*x1ry1r5.*y1r2-x1r4.*x2r2.*x1ry1r1.*y1r5-x1r4.*x2r5.*x1ry1r2.*y1r1-x1r5.*x2r1.*x1ry1r2.*y1r4-x1r5.*x2r2.*x1ry1r4.*y1r1-x1r5.*x2r4.*x1ry1r1.*y1r2-x1r1.*x2r3.*x1ry1r4.*y1r5-x1r1.*x2r4.*x1ry1r5.*y1r3-x1r1.*x2r5.*x1ry1r3.*y1r4-x1r3.*x2r1.*x1ry1r5.*y1r4-x1r3.*x2r4.*x1ry1r1.*y1r5-x1r3.*x2r5.*x1ry1r4.*y1r1-x1r4.*x2r1.*x1ry1r3.*y1r5-x1r4.*x2r3.*x1ry1r5.*y1r1-x1r4.*x2r5.*x1ry1r1.*y1r3-x1r5.*x2r1.*x1ry1r4.*y1r3-x1r5.*x2r3.*x1ry1r1.*y1r4-x1r5.*x2r4.*x1ry1r3.*y1r1-x1r2.*x2r3.*x1ry1r5.*y1r4-x1r2.*x2r4.*x1ry1r3.*y1r5-x1r2.*x2r5.*x1ry1r4.*y1r3-x1r3.*x2r2.*x1ry1r4.*y1r5-x1r3.*x2r4.*x1ry1r5.*y1r2-x1r3.*x2r5.*x1ry1r2.*y1r4-x1r4.*x2r2.*x1ry1r5.*y1r3-x1r4.*x2r3.*x1ry1r2.*y1r5-x1r4.*x2r5.*x1ry1r3.*y1r2-x1r5.*x2r2.*x1ry1r3.*y1r4-x1r5.*x2r3.*x1ry1r4.*y1r2-x1r5.*x2r4.*x1ry1r2.*y1r3).*-2.0,t75.*(t304+t305+t306+t307+t308+t309+t310+t311+t312+t313+t314+t315-x1r2.*x2r3.*x1ry1r5.*y1r4-x1r2.*x2r4.*x1ry1r3.*y1r5-x1r2.*x2r5.*x1ry1r4.*y1r3-x1r3.*x2r2.*x1ry1r4.*y1r5-x1r3.*x2r4.*x1ry1r5.*y1r2-x1r3.*x2r5.*x1ry1r2.*y1r4-x1r4.*x2r2.*x1ry1r5.*y1r3-x1r4.*x2r3.*x1ry1r2.*y1r5-x1r4.*x2r5.*x1ry1r3.*y1r2-x1r5.*x2r2.*x1ry1r3.*y1r4-x1r5.*x2r3.*x1ry1r4.*y1r2-x1r5.*x2r4.*x1ry1r2.*y1r3).*2.0,t75.*(t316+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327-x1r1.*x2r3.*x1ry1r4.*y1r5-x1r1.*x2r4.*x1ry1r5.*y1r3-x1r1.*x2r5.*x1ry1r3.*y1r4-x1r3.*x2r1.*x1ry1r5.*y1r4-x1r3.*x2r4.*x1ry1r1.*y1r5-x1r3.*x2r5.*x1ry1r4.*y1r1-x1r4.*x2r1.*x1ry1r3.*y1r5-x1r4.*x2r3.*x1ry1r5.*y1r1-x1r4.*x2r5.*x1ry1r1.*y1r3-x1r5.*x2r1.*x1ry1r4.*y1r3-x1r5.*x2r3.*x1ry1r1.*y1r4-x1r5.*x2r4.*x1ry1r3.*y1r1).*2.0,t75.*(t328+t329+t330+t331+t332+t333+t334+t335+t336+t337+t338+t339-x1r1.*x2r2.*x1ry1r5.*y1r4-x1r1.*x2r4.*x1ry1r2.*y1r5-x1r1.*x2r5.*x1ry1r4.*y1r2-x1r2.*x2r1.*x1ry1r4.*y1r5-x1r2.*x2r4.*x1ry1r5.*y1r1-x1r2.*x2r5.*x1ry1r1.*y1r4-x1r4.*x2r1.*x1ry1r5.*y1r2-x1r4.*x2r2.*x1ry1r1.*y1r5-x1r4.*x2r5.*x1ry1r2.*y1r1-x1r5.*x2r1.*x1ry1r2.*y1r4-x1r5.*x2r2.*x1ry1r4.*y1r1-x1r5.*x2r4.*x1ry1r1.*y1r2).*2.0,t75.*(t340+t341+t342+t343+t344+t345+t346+t347+t348+t349+t350+t351-x1r1.*x2r2.*x1ry1r3.*y1r5-x1r1.*x2r3.*x1ry1r5.*y1r2-x1r1.*x2r5.*x1ry1r2.*y1r3-x1r2.*x2r1.*x1ry1r5.*y1r3-x1r2.*x2r3.*x1ry1r1.*y1r5-x1r2.*x2r5.*x1ry1r3.*y1r1-x1r3.*x2r1.*x1ry1r2.*y1r5-x1r3.*x2r2.*x1ry1r5.*y1r1-x1r3.*x2r5.*x1ry1r1.*y1r2-x1r5.*x2r1.*x1ry1r3.*y1r2-x1r5.*x2r2.*x1ry1r1.*y1r3-x1r5.*x2r3.*x1ry1r2.*y1r1).*2.0,t75.*(t352+t353+t354+t355+t356+t357+t358+t359+t360+t361+t362+t363-x1r1.*x2r2.*x1ry1r4.*y1r3-x1r1.*x2r3.*x1ry1r2.*y1r4-x1r1.*x2r4.*x1ry1r3.*y1r2-x1r2.*x2r1.*x1ry1r3.*y1r4-x1r2.*x2r3.*x1ry1r4.*y1r1-x1r2.*x2r4.*x1ry1r1.*y1r3-x1r3.*x2r1.*x1ry1r4.*y1r2-x1r3.*x2r2.*x1ry1r1.*y1r4-x1r3.*x2r4.*x1ry1r2.*y1r1-x1r4.*x2r1.*x1ry1r2.*y1r3-x1r4.*x2r2.*x1ry1r3.*y1r1-x1r4.*x2r3.*x1ry1r1.*y1r2).*2.0];
end