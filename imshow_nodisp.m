function varargout = imshow_nodisp(im)
% An IMSHOW implementation that works even when Matlab runs -nodisplay.
%
% from http://stackoverflow.com/a/28717260/1114328
%
% Only we don't scale the figure window to reflect the image size. Consequently
% the ugly pixel interpolation is directly apparent. IMSHOW has it too, but it
% tries to hide it by scaling the figure window at once.
%
% Input arguments:
%  IM  HxWxD image.
%
% Output arguments:
%  HND  Handle to the drawn image (optional).
%
  [h,w,~] = size(im);

  x = [0 w; 0 w] + 0.5;
  y = [0 0; h h] + 0.5;
  z = [0 0; 0 0];

  hnd = surf(x, y, z, flipud(im), 'FaceColor','texturemap', 'EdgeColor','none');

  view(2);
  axis equal tight off;

  if nargout > 0
    varargout = hnd;
  end
end

