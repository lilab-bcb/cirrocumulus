import {format} from 'd3-format';

export const numberFormat = format('.1f');
export const numberFormat2f = format('.2f');

export const numberFormat2g = format('.2g');
export const intFormat = format(',');
export const numberFormat0 = format('.0f');

export function formatNumber(value) {
  let s = numberFormat2f(value);
  if (s.endsWith('.00')) {
    s = s.substring(0, s.lastIndexOf('.'));
    if (s === '0' && value !== 0) {
      s = numberFormat2g(value);
    }
  }
  return s;
}
