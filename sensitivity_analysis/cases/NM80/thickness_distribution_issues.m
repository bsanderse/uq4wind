% thickness issues in NM80 case

Z = [
0.00
2.00
4.00
6.00
8.00
10.00
12.00
14.00
16.00
18.00
20.00
22.00
24.00
26.00
28.00
30.00
32.00
34.00
36.00
37.00
38.00
38.40
38.80];

TC = [
99.99
96.41
80.53
65.08
51.67
40.30
32.53
28.40
25.62
23.77
22.25
20.99
20.03
19.40
19.03
18.79
18.60
18.39
17.95
17.39
16.33
15.70
14.84];

Z_airfoil = [
    11.876
    17.820
    28.976
    35.535];

TC_airfoil = [
    33.3
    24.3
    19.7
    18.7];

figure
plot(Z,TC,'x-','LineWidth',2,'markersize',10)
hold on
h2 = plot(Z_airfoil,TC_airfoil,'o','markersize',10);
set(h2, 'markerfacecolor', get(h2, 'color')); % Use same color to fill in markers


xlim([0 40])
ylim([0 50])
xlabel('r [m]');
ylabel('t/c [%]');
grid
legend('planform data','airfoil sections')
set(gca,'FontSize',14)

