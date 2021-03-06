if (MPC_TrajectoryType == 0)
    [MPC_Trajectory, MPC_VisualizationLimits] = GenerateTestTrajectory;
elseif (MPC_TrajectoryType == 1)
    [MPC_Trajectory, MPC_VisualizationLimits] = GenerateTestTrajectory_FigureEight;
end

% Simulated obstacle (defined in inertial frame)
ObstacleAvoidanceEnabled = true;
Obstacle_Avoidance_Clearance = 0.05;
RandomizedObstacles = 0; % enable random obstacles by setting this >0

MPC_Obstacles = [];
if (MPC_EnableStaticObstacles)    
    % Define some static obstacles that are always present
    %             X     Y     R   
    MPC_Obstacles = [ 4,   1.9,  0.2;
                 5.5,  1.6,  0.8;
                 3.10, 2.2,  0.3;
                 5.7, -1.93, 1.0;
                 6.6,  0.3,  0.4;
                -5.3,  1.8,  0.4;
                -6.1, -1.0, 0.6];
end
            
if (MPC_RandomObstacles > 0)
    for (n = 1:MPC_RandomObstacles)
        obs_x = MPC_VisualizationLimits.x_min + (MPC_VisualizationLimits.x_max-MPC_VisualizationLimits.x_min)*rand(1,1);
        obs_y = MPC_VisualizationLimits.y_min + (MPC_VisualizationLimits.y_max-MPC_VisualizationLimits.y_min)*rand(1,1);
        obs_r = 0.2 + 0.8 * rand(1,1); % between 0.2-1.0 meter radius
        MPC_Obstacles = [MPC_Obstacles;
                        [obs_x, obs_y, obs_r]];
    end
end

% If obstacle avoidance should not be tested, remove all obstacles
% Conversely add null-obstacles since at least 5 obstacles should be in the obstacle list
if (MPC_EnableStaticObstacles == false && MPC_RandomObstacles < 5)
    MPC_Obstacles = [MPC_Obstacles; repmat([100, 100, 1], [(5-size(MPC_Obstacles,1)),1])]; % needs at least 5 obstacles, here just set to far away
end  