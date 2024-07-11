using PyPlot

# 设置后端为 Qt5Agg 以确保弹出窗口
PyPlot.ion()  # 打开交互模式

# 绘制图像
x = 1:10
y = rand(10)
plot(x, y)
xlabel("X Axis")
ylabel("Y Axis")
title("Sample Plot")

# 显示图像
PyPlot.show(block=true)

# PyPlot.pause(100.0)