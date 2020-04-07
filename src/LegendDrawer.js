export function getCategoricalLegendSize(context, name, categories) {
    context.font = '14px Helvetica';
    let maxWidth = context.measureText(name).width;
    categories.forEach(value => maxWidth = Math.max(maxWidth, context.measureText(value).width));
    return {width: maxWidth + 14, height: categories.length * 12 + 4};
}

export function drawCategoricalLegend(context, scale, name, categories) {

    context.font = '14px Helvetica';
    context.textAlign = 'left';
    context.textBaseline = 'bottom';
    const height = 12;
    for (let i = 0; i < categories.length; i++) {
        const category = categories[i];
        context.fillStyle = scale(category);
        context.fillRect(0, i * height, 10, 10);
        context.fillStyle = 'black';
        context.fillText('' + category, 12, i * height + height);
    }
}


