import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_doubly_reinforced_beam(height, width, cover, topDiameter, topCount, botDiameter, botCount, mu, ax=None):
    """
    Draw a doubly reinforced concrete beam cross‐section with randomly oriented fibers.
    
    Parameters:
      height      : Overall height of the beam.
      width       : Overall width of the beam.
      cover       : Concrete cover to the reinforcement.
      topDiameter : Diameter of top rebars.
      topCount    : Number of top rebars.
      botDiameter : Diameter of bottom rebars.
      botCount    : Number of bottom rebars.
      mu          : A factor used to determine the number of fibers.
      ax          : (Optional) Matplotlib Axes object to draw into.
                    If not provided, a new figure and axes are created.
    """
    # If no axis is provided, create one and show the figure at the end.
    if ax is None:
        fig, ax = plt.subplots()
        show_fig = True
    else:
        show_fig = False

    ax.clear()
    ax.set_aspect('equal')
    
    # Define colors
    concreteColor = (0.7, 0.7, 0.7)  # Gray
    rebarColor    = (0.8, 0.2, 0.2)  # Reddish
    fiberColor    = (0, 0, 0)         # Black
    
    # Draw concrete cross section
    concrete_rect = patches.Rectangle((0, 0), width, height,
                                      facecolor=concreteColor, edgecolor='k')
    ax.add_patch(concrete_rect)
    
    # Calculate rebar center positions
    topRebarY = height - cover
    botRebarY = cover

    # Draw random fibers if applicable.
    # Determine number and length of fibers.
    if width > 74 or height > 74:
        numFibers = int(mu * 300)
        fiberLength = 10  # e.g., 10 mm (approx. 1 inch)
    else:
        numFibers = int(mu * (height * width) * 1.1)
        fiberLength = 1   # e.g., 1 unit
    
    fiberWidth = 0.2  # Thin fibers
    
    for j in range(numFibers):
        # Generate a random angle in radians (0 to 2*pi)
        angle = np.random.rand() * 2 * np.pi
        deltaX = np.cos(angle) * fiberLength
        deltaY = np.sin(angle) * fiberLength

        # Compute allowable start coordinates so that the entire fiber lies within the beam
        xMin = max(0, -deltaX)
        yMin = max(0, -deltaY)
        xMax = min(width, width - deltaX)
        yMax = min(height, height - deltaY)

        xStart = xMin + (xMax - xMin) * np.random.rand()
        yStart = yMin + (yMax - yMin) * np.random.rand()

        xEnd = xStart + deltaX
        yEnd = yStart + deltaY

        ax.plot([xStart, xEnd], [yStart, yEnd], color=fiberColor, linewidth=fiberWidth)
    
    # Draw top rebars
    if topCount > 1:
        topSpacing = (width - 2 * (cover + topDiameter / 2)) / (topCount - 1)
        for i in range(topCount):
            x = cover + topDiameter / 2 + i * topSpacing
            circle = patches.Circle((x, topRebarY), radius=topDiameter / 2,
                                    facecolor=rebarColor, edgecolor='k')
            ax.add_patch(circle)
    elif topCount == 1:
        x = width / 2
        circle = patches.Circle((x, topRebarY), radius=topDiameter / 2,
                                facecolor=rebarColor, edgecolor='k')
        ax.add_patch(circle)
    
    # Draw bottom rebars
    if botCount > 1:
        botSpacing = (width - 2 * (cover + botDiameter / 2)) / (botCount - 1)
        for i in range(botCount):
            x = cover + botDiameter / 2 + i * botSpacing
            circle = patches.Circle((x, botRebarY), radius=botDiameter / 2,
                                    facecolor=rebarColor, edgecolor='k')
            ax.add_patch(circle)
    elif botCount == 1:
        x = width / 2
        circle = patches.Circle((x, botRebarY), radius=botDiameter / 2,
                                facecolor=rebarColor, edgecolor='k')
        ax.add_patch(circle)
    
    # Set plot limits with some extra room equal to cover
    ax.set_xlim([-cover, width + cover])
    ax.set_ylim([-cover, height + cover])
    
    ax.set_xlabel('Width')
    ax.set_ylabel('Height')
    ax.set_title('Beam Cross-Section')
    
    if show_fig:
        plt.show()
