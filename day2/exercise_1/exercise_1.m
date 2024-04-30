v_1 = [-1 0.4]';
v_2 = [1 0.7]';
v_3 = [0.5 -0.8]';
v = [v_1 v_2 v_3];

e = gram_schmid(v);

figure
hold on
grid on
axis equal
plot([0 v(1,1)], [0 v(2,1)], "Color","blue", "LineWidth",1)
plot([0 v(1,2)], [0 v(2,2)], "Color","blue", "LineWidth",1)
plot([0 v(1,3)], [0 v(2,3)], "Color","blue", "LineWidth",1)
plot([0 e(1,1)], [0 e(2,1)], "Color","red", "LineWidth",1)
plot([0 e(1,2)], [0 e(2,2)], "Color","red", "LineWidth",1)
xlim([-max(max(abs(v)))-0.5, max(max(abs(v)))+0.5])

v_1 = [-1 0.4 -2]';
v_2 = [1 0.7 0.3]';
v_3 = [0.5 -0.8 0.3]';
v = [v_1 v_2 v_3];

e = gram_schmid(v);

figure
hold on
grid on
axis equal
quiver3(0,0,0,v(1,1),v(2,1),v(3,1), "Color","blue", "LineWidth",1)
quiver3(0,0,0,v(1,2),v(2,2),v(3,2), "Color","blue", "LineWidth",1)
quiver3(0,0,0,v(1,3),v(2,3),v(3,3), "Color","blue", "LineWidth",1)
quiver3(0,0,0,e(1,1),e(2,1),e(3,1), "Color","red", "LineWidth",1)
quiver3(0,0,0,e(1,2),e(2,2),e(3,2), "Color","red", "LineWidth",1)
quiver3(0,0,0,e(1,3),e(2,3),e(3,3), "Color","red", "LineWidth",1)
xlim([-max(max(abs(v)))-0.5, max(max(abs(v)))+0.5])
