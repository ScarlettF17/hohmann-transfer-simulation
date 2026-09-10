import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Circle
from pathlib import Path
import time

"FUNCTION DEFINITIONS"

def accel(r, GM):
    Rsqr = np.sum(r**2)
    return -r * GM / Rsqr**(3/2)

def streak_full(frame, rMe, rV, rE, rM, rS, trail_len):
    # Update all object positions for the animation.
    start = max(0, frame - trail_len)

    pathMe.set_data(rMe[start:frame+1, 0], rMe[start:frame+1, 1])
    frontMe.set_data([rMe[frame, 0]], [rMe[frame, 1]])

    pathV.set_data(rV[start:frame+1, 0], rV[start:frame+1, 1])
    frontV.set_data([rV[frame, 0]], [rV[frame, 1]])

    pathE.set_data(rE[start:frame+1, 0], rE[start:frame+1, 1])
    frontE.set_data([rE[frame, 0]], [rE[frame, 1]])

    pathM.set_data(rM[start:frame+1, 0], rM[start:frame+1, 1])
    frontM.set_data([rM[frame, 0]], [rM[frame, 1]])

    pathS.set_data(rS[start:frame+1, 0], rS[start:frame+1, 1])
    frontS.set_data([rS[frame, 0]], [rS[frame, 1]])

    return pathMe, frontMe, pathV, frontV, pathE, frontE, pathM, frontM, pathS, frontS

"OUTPUT FOLDERS"

# Save files into clean GitHub-friendly folders.
project_folder = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
figures_folder = project_folder / "figures"
animation_folder = project_folder / "animation"

figures_folder.mkdir(exist_ok=True)
animation_folder.mkdir(exist_ok=True)

"VARIABLES AND CONSTANTS"

M_Sun = 1.9891e30
G = 1.4881e-34
GM = G * M_Sun

dt = 0.1
T = np.arange(0.0, 1825.0, dt)   # 5 years
N = T.size

"FOUND PLANET DATA"

# Perihelion orbital data (AU and AU/day).
rMerc_peri = 0.3075
vMerc_peri = 0.0340

rVen_peri = 0.7184
vVen_peri = 0.0202

rEarth_peri = 0.9830
vEarth_peri = 0.0175

rMars_peri = 1.3810
vMars_peri = 0.0153

# Arrays for position and velocity.
rMe = np.zeros((N, 2))
vMe = np.zeros_like(rMe)

rV = np.zeros((N, 2))
vV = np.zeros_like(rV)

rE = np.zeros((N, 2))
vE = np.zeros_like(rE)

rM = np.zeros((N, 2))
vM = np.zeros_like(rM)

rS = np.zeros((N, 2))
vS = np.zeros_like(rS)

"INITIAL CONDITIONS"

# Mercury at 9 o'clock.
rMe[0, :] = [-rMerc_peri, 0.0]
vMe[0, :] = [0.0, -vMerc_peri]

# Venus at 12 o'clock.
rV[0, :] = [0.0, rVen_peri]
vV[0, :] = [-vVen_peri, 0.0]

# Earth at 3 o'clock.
rE[0, :] = [rEarth_peri, 0.0]
vE[0, :] = [0.0, vEarth_peri]

# Mars at 6 o'clock.
rM[0, :] = [0.0, -rMars_peri]
vM[0, :] = [vMars_peri, 0.0]

# Satellite starts with Earth.
rS[0, :] = rE[0, :]
vS[0, :] = vE[0, :]

# Hohmann transfer parameters.
launch_day = 460.7
launch_factor = 1.072

landing_day = 695.3
landing_factor = 1.146

launch_index = int(launch_day / dt)
landing_index = int(landing_day / dt)

"EULER SIMULATION LOOP"

clockStart = time.time()

for t in range(N - 1):
    # Mercury
    aMe = accel(rMe[t, :], GM)
    vMe[t+1, :] = vMe[t, :] + aMe * dt
    rMe[t+1, :] = rMe[t, :] + vMe[t+1, :] * dt

    # Venus
    aV = accel(rV[t, :], GM)
    vV[t+1, :] = vV[t, :] + aV * dt
    rV[t+1, :] = rV[t, :] + vV[t+1, :] * dt

    # Earth
    aE = accel(rE[t, :], GM)
    vE[t+1, :] = vE[t, :] + aE * dt
    rE[t+1, :] = rE[t, :] + vE[t+1, :] * dt

    # Mars
    aM = accel(rM[t, :], GM)
    vM[t+1, :] = vM[t, :] + aM * dt
    rM[t+1, :] = rM[t, :] + vM[t+1, :] * dt

    # Satellite
    aS = accel(rS[t, :], GM)
    vS[t+1, :] = vS[t, :] + aS * dt

    # First burn: launch
    if t == launch_index:
        vS[t+1, :] = launch_factor * vS[t+1, :]

    # Second burn: arrival adjustment
    if t == landing_index:
        vS[t+1, :] = landing_factor * vS[t+1, :]

    rS[t+1, :] = rS[t, :] + vS[t+1, :] * dt

print("Duration of for loop: %5f s" % (time.time() - clockStart))

"B1 CALCULATIONS"

# Distance between satellite and Mars at every time step.
sep_SM = np.sqrt(np.sum((rS - rM)**2, axis=1))

# Separation on the selected landing day.
landing_sep = sep_SM[landing_index]

# Closest approach.
closest_index = np.argmin(sep_SM)
closest_day = T[closest_index]
closest_sep = sep_SM[closest_index]

# Distance of the satellite from the Sun at every time step.
rS_mag = np.sqrt(np.sum(rS**2, axis=1))

# First aphelion after launch.
aphelion_index = launch_index
for i in range(launch_index + 1, N - 1):
    if rS_mag[i - 1] < rS_mag[i] and rS_mag[i] > rS_mag[i + 1]:
        aphelion_index = i
        break

aphelion_day = T[aphelion_index]
aphelion_radius = rS_mag[aphelion_index]

print("\nB1 Results")
print("Launch Day = %.1f days" % launch_day)
print("Launch Speed Factor = %.3f" % launch_factor)
print("Landing Day = %.1f days" % landing_day)
print("Landing Speed Factor = %.3f" % landing_factor)
print("Separation on Landing Day = %.6e AU" % landing_sep)
print("Closest approach = %.6e AU on day %.1f" % (closest_sep, closest_day))
print("Satellite aphelion = %.6f AU on day %.1f" % (aphelion_radius, aphelion_day))

"PLOT: COMPLETED TRANSFER"

fig_save = plt.figure(figsize=(12, 12))
ax_save = fig_save.add_subplot(1, 1, 1)

ax_save.set_aspect("equal")
ax_save.set_xlim(-2.0, 2.0)
ax_save.set_ylim(-2.0, 2.0)

# Draw the Sun.
ax_save.add_artist(Circle((0.0, 0.0), 0.07, color="#ffaa00"))

# Plot planetary orbits.
ax_save.plot(rMe[:, 0], rMe[:, 1], label="Mercury")
ax_save.plot(rV[:, 0], rV[:, 1], label="Venus")
ax_save.plot(rE[:, 0], rE[:, 1], label="Earth")
ax_save.plot(rM[:, 0], rM[:, 1], label="Mars")

# Plot the satellite transfer orbit.
ax_save.plot(rS[:, 0], rS[:, 1], "m--", linewidth=2, label="Satellite")

# Mark launch and landing positions.
ax_save.plot(rS[launch_index, 0], rS[launch_index, 1], "mo", label="Launch")
ax_save.plot(rS[landing_index, 0], rS[landing_index, 1], "ko", label="Landing")
ax_save.plot(rM[landing_index, 0], rM[landing_index, 1], "ro")

details = (
    "Hohmann Transfer Details\n"
    f"Launch Day = {launch_day:.1f} days\n"
    f"Launch Speed Factor = {launch_factor:.3f}\n"
    f"Landing Day = {landing_day:.1f} days\n"
    f"Landing Speed Factor = {landing_factor:.3f}\n"
    f"Separation on Landing Day = {landing_sep:.6e} AU"
)

ax_save.text(-1.95, 1.75, details, fontsize=10, bbox=dict(facecolor="white", edgecolor="black", alpha=0.85))
ax_save.set_title("Earth to Mars Hohmann Transfer")
ax_save.set_xlabel("x (AU)")
ax_save.set_ylabel("y (AU)")
ax_save.legend(loc="lower left")

plt.tight_layout()
plt.savefig(figures_folder / "transfer.png", dpi=300)

"ANIMATION DATA"

# Use every 10th point so the animation runs faster.
rMe_anim = rMe[::10, :]
rV_anim = rV[::10, :]
rE_anim = rE[::10, :]
rM_anim = rM[::10, :]
rS_anim = rS[::10, :]

"ANIMATION"

fig1 = plt.figure(figsize=(12, 12))
ax1 = fig1.add_subplot(1, 1, 1)

ax1.set_aspect("equal")
ax1.set_xlim(-2.0, 2.0)
ax1.set_ylim(-2.0, 2.0)
ax1.set_title("Animation: Full System")

# Draw the Sun.
ax1.add_artist(Circle((0.0, 0.0), 0.1, color="#ffaa00"))

# Empty artists updated frame by frame.
pathMe, = ax1.plot([], [], "-", linewidth=1)
frontMe, = ax1.plot([], [], "o", markersize=5, label="Mercury")

pathV, = ax1.plot([], [], "-", linewidth=1)
frontV, = ax1.plot([], [], "o", markersize=6, label="Venus")

pathE, = ax1.plot([], [], "-", linewidth=1)
frontE, = ax1.plot([], [], "bo", markersize=7, label="Earth")

pathM, = ax1.plot([], [], "-", linewidth=1)
frontM, = ax1.plot([], [], "ro", markersize=7, label="Mars")

pathS, = ax1.plot([], [], "m--", linewidth=1.5)
frontS, = ax1.plot([], [], "m^", markersize=6, label="Satellite")

ax1.legend(loc="lower left")

ani1 = animation.FuncAnimation(
    fig1,
    streak_full,
    frames=len(rE_anim),
    repeat=False,
    fargs=(rMe_anim, rV_anim, rE_anim, rM_anim, rS_anim, 300),
    interval=20,
    blit=True
)

"SAVE ANIMATION"

print("\nSaving GIF... this may take a little time.")
ani1.save(animation_folder / "hohmann_transfer.gif", writer="pillow", fps=30)
print("GIF saved to:", animation_folder / "hohmann_transfer.gif")

plt.show()