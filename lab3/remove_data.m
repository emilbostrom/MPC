function remove_data(col)
  figure(1)
  h=findobj(gcf,'color',col);
  delete(h)
  figure(2)
  h=findobj(gcf,'color',col);
  delete(h)