# GeoDraw.jl
The idea for this project started nearly a year ago, when I had to make some illustrations related to satellite geodesy (a classical diagram of orbital elements). The standard procedure is to use TikZ. At that time, I could not, I was very new to using TikZ, and it turned out very daunting, especially the annotation aspect. I explored the Python route at that time as well, but gave up and moved on. Recently at EGU26, I saw a poster from the team StrataPy (https://github.com/Jack6228/stratapy) and was really inspired to revisit my problem of generating programmatic illustrations.

My Julia skills have improved over the last year; hopefully, I can make a working prototype.

The idea has been to make "primitives" that can be composed to create an illustration. The idea is to realise a 3D version with Makie, which one can fiddle with by changing the camera view and selecting an appropriate camera angle to export a .svg graphic. The geometry, style, composition, and annotations should be handled separately to enhance modularity and extensibility.

The idea would be to provide a contextual, easy-to-use abstraction for TikZ-like graphic primitives in Julia. So instead of a group of 3 arrows, it can be abstracted to the geodesy context of a frame with parameters such as origin, orientation, and scale. This also separates the concern of scientific context from the aspect of making it appear on the screen. The hope would be to get a WYSWIG, i.e., what you Science - what you get kind of system.

Let's see how much I can progress. Open to comments, suggestions or criticism of the idea or implementation.