function data = icra17_process(mode, data, mass, H_vic2bod, H_m402bod, H_bal2imu)

    materials = data.keys;

    switch mode
        case 'vicon'
            data = vicon(data, materials, mass, H_vic2bod, H_m402bod, H_bal2imu);
            
        case 'bluefox'
            data = bluefox(data, materials, mass, H_vic2bod, H_m402bod, H_bal2imu);
            
        case 'both'
            data = vicon(data, materials, mass, H_vic2bod, H_m402bod, H_bal2imu);
            data = bluefox(data, materials, mass, H_vic2bod, H_m402bod, H_bal2imu);
    end
end

function data = vicon(data, materials, mass, H_vic2bod, H_m402bod, H_bal2imu)
    for m = 1:length(materials)
        ep = data(materials{m});

        [~,~,~,~,~,~, vei, ai, ~,~,~, iws] = process_stick(ep.v, ep.int, ep.acc, mass, [0;0;0], H_vic2bod, H_m402bod, H_bal2imu, -ep.off);
        ep.vei = vei;
        ep.ai = ai;
        ep.iws = iws;

        data(materials{m}) = ep;
    end
end

function data = bluefox(data, materials, mass, H_vic2bod, H_m402bod, H_bal2imu)
    for m = 1:length(materials)
        ep = data(materials{m});

        [~,~,~,~,~,~, vei, ai, ~,~,~, iws] = process_stick(ep.motrak, ep.int, ep.acc, mass, [0;0;0], H_vic2bod, H_m402bod, H_bal2imu);
        ep.bvei = vei;
        ep.bai = ai;
        ep.biws = iws;

        data(materials{m}) = ep;
    end
end
